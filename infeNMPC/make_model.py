"""
Pyomo model construction for infeNMPC.

Provides functions to build steady-state, finite-horizon, and
infinite-horizon NMPC models.  Model equations and variable declarations
are loaded via ``_get_model(options)``, which resolves ``options.model_module``
directly (primary) or falls back to ``_load_model`` by name (legacy).

Non-collocation cleanup is handled by passing ``clean_model='delete'`` to
``dae.collocation.apply_to`` (requires ``oag00002/pyomo`` branch ``mpc-rewrite``).
``DynamicModelInterface`` is created with ``clean_model=True`` so that all
MPC data-access operations use the safe sparse-access paths from
``contrib.mpc.interfaces.*_clean``.
"""
from .infNMPC_options import _import_settings  # noqa: F401  (used by legacy callers)
import pyomo.environ as pyo
from .model_equations import _get_model
from pyomo.contrib.mpc import DynamicModelInterface, ScalarData
from pyomo.dae import ContinuousSet, DerivativeVar
from pyomo.common.collections import ComponentMap
from pyomo.core.expr.visitor import replace_expressions
from pyomo.core.expr.visitor import identify_variables
from pyomo.opt import TerminationCondition
from .tools.indexing_tools import (
    _get_variable_key_for_data,
    _add_time_indexed_expression,
    _get_derivative_and_state_vars,
)
from .tools.collocation_tools import (
    _compute_gamma_from_collocation,
    _get_last_fe_pts,
    _lagrange_coeffs_at_endpoint,
)
from .tools.initialization_tools import _check_ic_consistency, _build_ic_data
from .tools.debug_tools import _report_constraint_violations


# ---------------------------------------------------------------------------
# Terminal-constraint helpers
# ---------------------------------------------------------------------------

def _terminal_vars(m, options):
    """
    Return an ordered list of variable names to include in the terminal
    constraint, based on ``options.terminal_constraint_variables``.

    Parameters
    ----------
    m : pyo.Block
        Block that has ``CV_index`` and ``MV_index``.
    options : Options

    Returns
    -------
    list of str
    """
    tvars = options.terminal_constraint_variables
    result = []
    if 'cv' in tvars:
        result += list(m.CV_index)
    if 'mv' in tvars:
        result += list(m.MV_index)
    return result


def _terminal_weights(options, var_names, stage_cost_index):
    """
    Return per-variable weights for the terminal constraint penalty.

    Weights come from ``options.stage_cost_weights``, indexed by each
    variable's position in *stage_cost_index* (the same ordering used for the
    stage cost).  Variables not found in *stage_cost_index* get weight 1.0.

    Parameters
    ----------
    options : Options
    var_names : list of str
    stage_cost_index : list of str
        Ordered list ``[cv0, cv1, ..., mv0, mv1, ...]`` matching the stage
        cost weight ordering.

    Returns
    -------
    list of float
    """
    c = options.stage_cost_weights or []
    pos = {v: i for i, v in enumerate(stage_cost_index)}
    return [c[pos[v]] if v in pos and pos[v] < len(c) else 1.0
            for v in var_names]


def _fix_initial_conditions_at_t0(block):
    """
    Fix initial conditions at t=0 for differential state variables.

    All equations_write constraint entries at t=0 are deleted before this is
    called (see ``_finite_block_gen``), so there are no algebraic coupling
    constraints at t=0.  Every differential state variable is fixed directly
    at its current value (loaded from initial_data by the caller).

    Parameters
    ----------
    block : pyo.ConcreteModel or pyo.Block
        Must have ``.time`` and ``.state_vars``.  Must have been discretised
        and had t=0 equations_write entries removed before this is called.
    """
    t0 = block.time.first()
    for var in block.state_vars:
        for index in var:
            time_val = index[-1] if isinstance(index, tuple) else index
            if time_val == t0:
                var[index].fix()


def _check_optimal(results, label=""):
    """Raise RuntimeError if the solver did not return an optimal solution."""
    tc = results.solver.termination_condition
    if tc != TerminationCondition.optimal:
        suffix = f" ({label})" if label else ""
        raise RuntimeError(
            f"Solver terminated non-optimally{suffix}: {tc}"
        )


def _ipopt_solver():
    """Return a pre-configured IPOPT solver instance for cold (initial) solves."""
    solver = pyo.SolverFactory('ipopt_v2')
    solver.options['linear_solver'] = 'ma57'
    solver.options['OF_ma57_automatic_scaling'] = 'yes'
    solver.options['max_iter'] = 6000
    solver.options['halt_on_ampl_error'] = 'yes'
    solver.options['bound_relax_factor'] = 0
    return solver


def _ipopt_warm_solver():
    """Return a pre-configured IPOPT solver instance for warm-started MPC updates.

    After a shift-and-load the primal/dual solution from the previous iteration
    is already a good starting point.  The warm-start settings keep IPOPT close
    to that point:

    * ``warm_start_init_point`` — tells IPOPT to accept the primal/dual values
      already in the model as the initial iterate rather than re-projecting them.
    * ``mu_init`` — sets the initial barrier parameter to a small value so that
      IPOPT starts in the neighbourhood of the primal-dual solution rather than
      from the default (usually 0.1).
    * ``bound_push`` / ``bound_frac`` — control how aggressively IPOPT pushes the
      starting point away from bounds.  Defaults (0.01) are appropriate for cold
      starts; much smaller values (1e-8) keep the warm-started iterate where it is.
    """
    solver = pyo.SolverFactory('ipopt_v2')
    solver.options['linear_solver'] = 'ma57'
    solver.options['OF_ma57_automatic_scaling'] = 'yes'
    solver.options['max_iter'] = 6000
    solver.options['halt_on_ampl_error'] = 'yes'
    solver.options['bound_relax_factor'] = 0
    # warm_start_init_point is intentionally omitted: Pyomo's ipopt_v2 interface
    # does not pass dual variable (multiplier) values to IPOPT without an explicit
    # Suffix setup.  Enabling it with zero duals gives IPOPT a badly-conditioned
    # initial KKT system and a misleading first search direction.  The primal warm
    # start is already fully effective: shift_values_by_time places the shifted
    # solution directly in the Pyomo model, and IPOPT reads those as its initial
    # primal iterate automatically.
    #
    # mu_init is intentionally omitted (IPOPT default = 0.1): the primal
    # infeasibility at the warm-started point grows as the Lyapunov constraint
    # becomes increasingly active with each iteration.  A fixed mu_init = 1e-3
    # is fine for nearly-feasible starts but prevents recovery when infeasibility
    # is large (e.g. 19.1 at iteration 2), causing IPOPT to spend all its
    # iterations making no progress.
    #
    # bound_push / bound_frac: tighter than IPOPT default (0.01) so that the
    # good warm-started variable values are preserved near their bounds rather
    # than being projected further into the interior at the start of each solve.
    solver.options['bound_push'] = 1e-8
    solver.options['bound_frac'] = 1e-8
    return solver


def _make_steady_state_model(m, options, fix_slacks=True):
    """
    Generate a steady-state model by fixing all derivatives to zero.

    Discretizes with a single finite element / collocation point so that the
    model represents a static NLP, then deactivates the discretization
    equations so the solver treats derivatives as parameters pinned to zero.

    Parameters
    ----------
    m : pyo.ConcreteModel
    options : Options
    fix_slacks : bool, optional
        If True (default) and the model declares ``m.slack_index``, fix all
        slack variables to 0 before solving.  Set to False for the IV solve
        so that the initial point is allowed to violate soft constraints.

    Returns
    -------
    pyo.ConcreteModel
    """
    _model = _get_model(options)

    m.time = ContinuousSet(bounds=(0, 1))
    m = _model.variables_initialize(m)

    if fix_slacks and hasattr(m, 'slack_index'):
        for sv_name in m.slack_index:
            sv_var = getattr(m, sv_name, None)
            if sv_var is not None and isinstance(sv_var, pyo.Var):
                for idx in sv_var.index_set():
                    sv_var[idx].fix(0)

    deriv_vars, state_vars = _get_derivative_and_state_vars(m)

    discretizer = pyo.TransformationFactory('dae.collocation')
    discretizer.apply_to(m, ncp=1, nfe=1, wrt=m.time, scheme='LAGRANGE-RADAU')

    m = _model.equations_write(m)

    for var in deriv_vars:
        for idx in var.index_set():
            time_val = idx[-1] if isinstance(idx, tuple) else idx
            if time_val in m.time:
                var[idx].fix(0)

    for var in deriv_vars:
        getattr(m, f"{var}_disc_eq").deactivate()

    return m


def _solve_steady_state_model(m, target, options, label="ss"):
    """
    Solve the steady-state model and return data at the initial time point.

    When *target* is ``None``, the objective is the average of the model's
    ``custom_objective`` stage cost at ``t=0`` and ``t=1`` (economic SS solve).
    Otherwise a standard penalty-from-target objective is used.

    Parameters
    ----------
    m : pyo.ConcreteModel
    target : ScalarData or None
    options : Options
    label : str
        Short name used as the stem of the pprint output file when
        ``options.model_output_dir`` is set
        (e.g. ``"ss"`` → ``<model_output_dir>/ss_model.txt``).

    Returns
    -------
    ScalarData
        Variable values at ``m.time.first()``.
    """
    ss_interface = DynamicModelInterface(m, m.time)

    if target is None:
        custom_objective = _get_model(options).custom_objective
        cost_fn = custom_objective(m, options)

        def ss_obj_rule(m):
            return (cost_fn(m, 0) + cost_fn(m, 1)) / 2

        m.objective = pyo.Objective(expr=ss_obj_rule(m), sense=pyo.minimize)

    else:
        var_set, tr_cost = ss_interface.get_penalty_from_target(target)
        m.target_set = var_set
        m.tracking_cost = tr_cost
        m.objective = pyo.Objective(expr=sum(m.tracking_cost[:, 0]))

    if options.model_output_dir:
        import os
        os.makedirs(options.model_output_dir, exist_ok=True)
        with open(os.path.join(options.model_output_dir, f"{label}_model.txt"), "w") as f:
            m.pprint(ostream=f)

    solver = _ipopt_solver()
    results = solver.solve(m, tee=options.tee_flag)
    try:
        _check_optimal(results, "steady-state")
    except RuntimeError:
        if options.debug_flag:
            _report_constraint_violations(m, label=f"{label} steady-state failure")
        raise
    if options.debug_flag:
        _report_constraint_violations(m, label=f"{label} steady-state")

    m.ss_obj_value = pyo.value(m.objective)

    return ss_interface.get_data_at_time(m.time.first())


def _make_infinite_horizon_model(m, options):
    """
    Construct and initialize a two-block infinite-horizon NMPC model.

    Builds a steady-state model to compute setpoint values, then constructs
    the finite and infinite horizon blocks and links them with interface
    constraints.

    Parameters
    ----------
    m : pyo.ConcreteModel
    options : Options

    Returns
    -------
    pyo.ConcreteModel
        Model with ``finite_block``, ``infinite_block``, and ``interface``
        attributes.
    """
    # --- Steady-state solve for setpoint values ---
    print('Writing Steady State Model')
    m_ss = pyo.ConcreteModel()
    m_ss = _make_steady_state_model(m_ss, options)

    if options.objective == 'economic' or options.tracking_setpoint == 'economic':
        m_ss_target = None
    else:
        m_ss_target = ScalarData({
            _get_variable_key_for_data(m_ss, var): pyo.value(m_ss.setpoints[var])
            for var in m_ss.setpoints.index_set()
        })

    print('Solving Steady State Model')
    steady_state_data = _solve_steady_state_model(m_ss, m_ss_target, options, label="ss")

    steady_state_values = {
        **{cv: steady_state_data.get_data_from_key(_get_variable_key_for_data(m_ss, cv))
           for cv in m_ss.CV_index},
        **{mv: steady_state_data.get_data_from_key(_get_variable_key_for_data(m_ss, mv))
           for mv in m_ss.MV_index},
    }

    # --- Initial condition data ---
    initial_data = _build_ic_data(options)

    print('Writing Infinite Horizon Model')

    m.finite_block = pyo.Block()
    m.infinite_block = pyo.Block()
    m.finite_block.steady_state_values = steady_state_values
    m.infinite_block.steady_state_values = steady_state_values

    if options.objective == 'economic':
        m.finite_block.ss_obj_value = m_ss.ss_obj_value
        m.infinite_block.ss_obj_value = m_ss.ss_obj_value

    m.finite_block = _finite_block_gen(m.finite_block, options)
    m.infinite_block = _infinite_block_gen(m.infinite_block, options)
    m = _link_blocks(m)

    # Lyapunov stability constraint for infinite-horizon models.
    # V = finite_block.phi_track (Riemann sum over FE endpoints)
    #   + infinite_block.phi_track[tau=1] (ODE integral over transformed domain)
    # Constraint: V - V_prev <= -delta * first_stage_cost_prev
    if options.lyap_flag:
        m.V_prev = pyo.Param(mutable=True, initialize=1e10)
        m.first_stage_cost_prev = pyo.Param(mutable=True, initialize=0.0)
        tau_last = m.infinite_block.time.last()
        lct = getattr(options, 'lyap_constraint_type', 'hard')
        if lct == 'hard':
            m.lyap_stability_constraint = pyo.Constraint(
                expr=(
                    m.finite_block.phi_track + m.infinite_block.phi_track[tau_last] - m.V_prev
                    <= -options.lyap_delta * m.first_stage_cost_prev
                )
            )
        elif lct == 'soft':
            m.lyap_slack = pyo.Var(initialize=0.0, domain=pyo.NonNegativeReals)
            m.lyap_stability_constraint = pyo.Constraint(
                expr=(
                    m.finite_block.phi_track + m.infinite_block.phi_track[tau_last] - m.V_prev
                    <= -options.lyap_delta * m.first_stage_cost_prev + m.lyap_slack
                )
            )
        # 'none': infrastructure built, no constraint added

    m.interface = DynamicModelInterface(
        m.finite_block, m.finite_block.time, clean_model=True
    )

    if options.initialize_with_initial_data:
        for time in m.finite_block.time:
            m.interface.load_data(initial_data, time_points=time)
    else:
        m.interface.load_data(initial_data, time_points=0)

    _fix_initial_conditions_at_t0(m.finite_block)
    _check_ic_consistency(m.finite_block)

    return m


def _make_finite_horizon_model(m, options):
    """
    Construct and initialize a finite-horizon NMPC model.

    Parameters
    ----------
    m : pyo.ConcreteModel
    options : Options

    Returns
    -------
    pyo.ConcreteModel
    """
    m_ss = pyo.ConcreteModel()
    m_ss = _make_steady_state_model(m_ss, options)

    if options.objective == 'economic' or options.tracking_setpoint == 'economic':
        m_ss_target = None
    else:
        m_ss_target = ScalarData({
            _get_variable_key_for_data(m_ss, var): pyo.value(m_ss.setpoints[var])
            for var in m_ss.setpoints.index_set()
        })

    steady_state_data = _solve_steady_state_model(m_ss, m_ss_target, options, label="ss")

    steady_state_values = {
        **{cv: steady_state_data.get_data_from_key(_get_variable_key_for_data(m_ss, cv))
           for cv in m_ss.CV_index},
        **{mv: steady_state_data.get_data_from_key(_get_variable_key_for_data(m_ss, mv))
           for mv in m_ss.MV_index},
    }

    initial_data = _build_ic_data(options)

    print('Writing Finite Horizon Model')

    if options.objective == 'economic':
        m.ss_obj_value = m_ss.ss_obj_value
    m.steady_state_values = steady_state_values
    m = _finite_block_gen(m, options)

    m.interface = DynamicModelInterface(m, m.time, clean_model=True)

    # Warm-start at t > t0 with steady-state values.
    # Skipping t=0 preserves the Pyomo initialize=... values there, which are
    # the model's intended initial conditions for state variables not covered by
    # initial_values (e.g. tray compositions / holdups in distillation models).
    # Interior and final time points still receive a physically consistent
    # steady-state warm-start, which is important for convergence on large models.
    t_first = m.time.first()
    for t in m.time:
        if t != t_first:
            m.interface.load_data(steady_state_data, time_points=t)

    if options.initialize_with_initial_data:
        for time in m.time:
            m.interface.load_data(initial_data, time_points=time)
    else:
        m.interface.load_data(initial_data, time_points=t_first)

    _fix_initial_conditions_at_t0(m)
    _check_ic_consistency(m)

    return m


def _finite_block_gen(m, options):
    """
    Populate a finite-horizon block with variables, discretization, and equations.

    Uses LAGRANGE-RADAU collocation with ``clean_model='delete'`` to remove
    spurious variable and constraint entries at non-collocation points.

    Parameters
    ----------
    m : pyo.ConcreteModel or pyo.Block
    options : Options

    Returns
    -------
    pyo.ConcreteModel or pyo.Block
    """
    _model = _get_model(options)

    m.time = ContinuousSet(bounds=(0, options.nfe_finite * options.sampling_time))
    m = _model.variables_initialize(m)

    deriv_vars, state_vars = _get_derivative_and_state_vars(m)
    m.deriv_vars = deriv_vars
    m.state_vars = state_vars

    discretizer = pyo.TransformationFactory('dae.collocation')
    discretizer.apply_to(
        m, ncp=options.ncp_finite, nfe=options.nfe_finite,
        wrt=m.time, scheme='LAGRANGE-RADAU', clean_model='delete',
    )
    if options.ncp_finite > 1:
        for mv in m.MV_index:
            m = discretizer.reduce_collocation_points(
                m, var=getattr(m, mv), ncp=1, contset=m.time
            )
        if hasattr(m, 'slack_index'):
            for sv in m.slack_index:
                if hasattr(m, sv):
                    m = discretizer.reduce_collocation_points(
                        m, var=getattr(m, sv), ncp=1, contset=m.time
                    )

    _model.equations_write(m)

    # RADAU: delete all constraint entries at t=0 written by equations_write.
    # t=0 is not a collocation point; algebraic variables (CVs, slacks, MVs)
    # have no ODE equation there.  equations_write must run post-discretization
    # (Pyomo DAE dictionary expansion only works after the time set is populated),
    # so it iterates over all of m.time including t=0.  Those t=0 entries must
    # not exist — they couple algebraic vars at the non-collocation starting
    # point and cause IC inconsistencies when state-var ICs are updated.
    _t0_del = m.time.first()
    for _con in m.component_objects(pyo.Constraint):
        _keys_at_t0 = [
            _idx for _idx in list(_con.keys())
            if (_idx[-1] if isinstance(_idx, tuple) else _idx) == _t0_del
        ]
        for _idx in _keys_at_t0:
            del _con[_idx]

    if options.lyap_flag:
        track_list_ft = list(m.CV_index) + list(m.MV_index)
        c_raw_ft = options.stage_cost_weights or []
        c_track_ft = (
            c_raw_ft[:len(track_list_ft)]
            if len(c_raw_ft) >= len(track_list_ft)
            else [1.0] * len(track_list_ft)
        )
        # Sum tracking cost at each FE right-endpoint (Radau points), skipping t=0.
        # With LAGRANGE-RADAU, FE endpoints are collocation points, so MVs exist there.
        fe_pts_ft = [t for t in m.time.get_finite_elements() if t != m.time.first()]

        def _phi_track_expr(m):
            return sum(
                c_track_ft[i] * (
                    _add_time_indexed_expression(m, var, t)
                    - m.steady_state_values[var]
                )**2
                for t in fe_pts_ft
                for i, var in enumerate(track_list_ft)
            )

        m.phi_track = pyo.Expression(rule=_phi_track_expr)

        if not options.infinite_horizon:
            # For finite-only: Lyapunov constraint lives here.
            # For infinite horizon: combined constraint is added on the parent ConcreteModel
            # in _make_infinite_horizon_model, referencing both blocks.
            m.V_prev = pyo.Param(mutable=True, initialize=1e10)
            m.first_stage_cost_prev = pyo.Param(mutable=True, initialize=0.0)
            lct = getattr(options, 'lyap_constraint_type', 'hard')
            if lct == 'hard':
                m.lyap_stability_constraint = pyo.Constraint(
                    expr=(
                        m.phi_track - m.V_prev
                        <= -options.lyap_delta * m.first_stage_cost_prev
                    )
                )
            elif lct == 'soft':
                m.lyap_slack = pyo.Var(initialize=0.0, domain=pyo.NonNegativeReals)
                m.lyap_stability_constraint = pyo.Constraint(
                    expr=(
                        m.phi_track - m.V_prev
                        <= -options.lyap_delta * m.first_stage_cost_prev + m.lyap_slack
                    )
                )
            # 'none': infrastructure built, no constraint added

    # --- Terminal constraints (finite-horizon controller only) ---
    # For infinite-horizon models, the terminal point is τ=1 of the infinite
    # block, NOT the end of the finite block, so we skip this here.
    if options.terminal_constraint_type in ('hard', 'soft') and not options.infinite_horizon:
        stage_cost_index = list(m.CV_index) + list(m.MV_index)
        vars_tc = _terminal_vars(m, options)
        weights = _terminal_weights(options, vars_tc, stage_cost_index)
        t_final = m.time.last()   # RADAU: right-endpoint IS a collocation point

        if options.terminal_constraint_type == 'hard':
            def _finite_tc_rule(m, var_name):
                return (
                    _add_time_indexed_expression(m, var_name, t_final)
                    == m.steady_state_values[var_name]
                )
            m.finite_terminal_constraints = pyo.Constraint(
                vars_tc, rule=_finite_tc_rule
            )

        else:  # 'soft'
            m.finite_terminal_soft_penalty = pyo.Expression(
                expr=sum(
                    w * (_add_time_indexed_expression(m, v, t_final)
                         - m.steady_state_values[v])**2
                    for v, w in zip(vars_tc, weights)
                )
            )

    return m


def _transform_model_derivatives(model, deriv_vars):
    """
    Replace ``dxdt[t]`` with ``gamma/(sampling_time) * (1 - t²) * dxdt[t]``
    in all user constraints of the infinite-horizon block.

    This maps the physical derivative to the compressed ``[0, 1]`` time
    domain used by the infinite-horizon ODE integration.  Discretization
    equations (``*_disc_eq``) and terminal cost constraints are left unchanged.

    Parameters
    ----------
    model : pyo.Block
        The infinite-horizon block (must have ``model.gamma`` and
        ``model.sampling_time`` attributes).
    deriv_vars : set
        Set of ``DerivativeVar`` components whose occurrences should be
        replaced.
    """
    gamma = model.gamma
    sampling_time = model.sampling_time

    for con in model.component_objects(pyo.Constraint, active=True, descend_into=True):
        if not con.is_indexed():
            continue
        if (con.name.endswith("_disc_eq") or con.name.endswith("terminal_cost")
                or con.name.endswith("phi_track_ode")):
            continue

        for index in list(con.keys()):
            expr = con[index].expr
            if expr is None:
                continue

            if isinstance(index, tuple):
                if len(index) == 0:
                    continue
                time_index = index[-1]
            else:
                time_index = index

            replacement_map = ComponentMap()
            for var in identify_variables(expr, include_fixed=False):
                if var.parent_component() in deriv_vars:
                    replacement_map[var] = (
                        gamma / sampling_time * (1 - time_index**2) * var
                    )

            if replacement_map:
                new_expr = replace_expressions(expr, substitution_map=replacement_map)
                con[index].set_value(new_expr)

    # NOTE: τ=1 entries are intentionally NOT cleaned up here.
    # phi_track[1] must always exist for the Lyapunov stability constraint,
    # and differential CVs need their τ=1 value for hard terminal constraints.


def _infinite_block_gen(m, options):
    """
    Build the infinite-horizon block of the NLP.

    Sets up the transformed time domain ``[0, 1]``, calls
    ``variables_initialize`` and ``equations_write``, appends the auxiliary
    ODE for the Lyapunov integral variable ``phi``, discretizes with
    LAGRANGE-LEGENDRE collocation and ``clean_model='delete'``, then applies
    the derivative substitution that maps physical time to ``[0, 1]``.

    Parameters
    ----------
    m : pyo.Block
    options : Options

    Returns
    -------
    pyo.Block
    """
    _model = _get_model(options)

    m.time = ContinuousSet(bounds=(0, 1))
    m = _model.variables_initialize(m)

    # --- Resolve gamma before building constraints that capture it ---
    # If gamma is not provided, choose it so that the first collocation
    # point tau_1 satisfies tau_1 = tanh(gamma * sampling_time), i.e.
    # gamma = atanh(tau_1) / sampling_time.
    if options.gamma is None:
        gamma = _compute_gamma_from_collocation(
            options.nfe_infinite, options.ncp_infinite, options.sampling_time
        )
        print(f"Auto-selected gamma = {gamma:.6g} "
              f"(tau_1 = first collocation point, sampling_time = {options.sampling_time})")
    else:
        gamma = options.gamma

    def _add_dphidt(m):
        # --- phi: primary terminal cost integral ---
        m.phi = pyo.Var(m.time, initialize=0, domain=pyo.NonNegativeReals)
        m.phi[0].fix(0)
        m.dphidt = DerivativeVar(m.phi, wrt=m.time)

        m.stage_cost_index = pyo.Set(
            initialize=lambda m: list(m.CV_index) + list(m.MV_index)
        )
        c = options.stage_cost_weights

        if options.objective == 'economic':
            custom_objective = _model.custom_objective
            cost_fn = custom_objective(m, options)

            def _custom_terminal_cost_rule(m, t):
                if t == 0 or t == 1:
                    return pyo.Constraint.Skip
                return (
                    (gamma / options.sampling_time * (1 - t**2))
                    * m.dphidt[t] == cost_fn(m, t) - m.ss_obj_value
                )

            m.terminal_cost = pyo.Constraint(m.time, rule=_custom_terminal_cost_rule)

        else:
            def _terminal_cost_rule(m, t):
                if t == 0 or t == 1:
                    return pyo.Constraint.Skip
                stage_cost = sum(
                    c[i] * (
                        _add_time_indexed_expression(m, var_name, t)
                        - m.steady_state_values[var_name]
                    )**2
                    for i, var_name in enumerate(m.stage_cost_index)
                )
                return (
                    (gamma / options.sampling_time * (1 - t**2))
                    * m.dphidt[t] == stage_cost
                )

            m.terminal_cost = pyo.Constraint(m.time, rule=_terminal_cost_rule)

        # --- phi_track: quadratic tracking Lyapunov integral over CVs and MVs ---
        track_list = list(m.CV_index) + list(m.MV_index)
        c_raw = options.stage_cost_weights or []
        c_track = c_raw[:len(track_list)] if len(c_raw) >= len(track_list) else [1.0] * len(track_list)

        m.phi_track = pyo.Var(m.time, initialize=0, domain=pyo.NonNegativeReals)
        m.phi_track[0].fix(0)
        m.dphidt_track = DerivativeVar(m.phi_track, wrt=m.time)

        def _phi_track_rule(m, t):
            if t == 0 or t == 1:
                return pyo.Constraint.Skip
            tracking_cost = sum(
                c_track[i] * (
                    _add_time_indexed_expression(m, var, t)
                    - m.steady_state_values[var]
                )**2
                for i, var in enumerate(track_list)
            )
            return (
                (gamma / options.sampling_time * (1 - t**2))
                * m.dphidt_track[t] == tracking_cost
            )

        m.phi_track_ode = pyo.Constraint(m.time, rule=_phi_track_rule)

        return m

    deriv_vars, state_vars = _get_derivative_and_state_vars(m)
    m = _add_dphidt(m)

    m.gamma = gamma
    m.sampling_time = options.sampling_time

    deriv_vars, state_vars = _get_derivative_and_state_vars(m)
    m.deriv_vars = deriv_vars
    m.state_vars = state_vars

    discretizer = pyo.TransformationFactory('dae.collocation')
    discretizer.apply_to(
        m, ncp=options.ncp_infinite, nfe=options.nfe_infinite,
        wrt=m.time, scheme='LAGRANGE-LEGENDRE', clean_model='delete',
    )

    _model.equations_write(m)

    # --- Terminal constraints for the infinite-horizon block ---
    #
    # LEGENDRE collocation places quadrature points strictly inside each
    # finite element — τ=1 (the global endpoint) is NOT a collocation point.
    # Pyomo's internal continuity equations DO pin τ=1 for differential state
    # variables (those with a DerivativeVar) by extrapolating the polynomial
    # from the last element.  Algebraic variables (CVs without a DerivativeVar
    # and all MVs) have no such continuity equation, so their τ=1 index is a
    # genuinely free variable.
    #
    # Strategy
    # --------
    # • Differential CVs  → direct access at τ=1 (already determined).
    # • Algebraic CVs/MVs → Lagrange extrapolation from the last FE's
    #                        collocation points.
    if options.terminal_constraint_type in ('hard', 'soft'):
        stage_cost_index = list(m.CV_index) + list(m.MV_index)
        vars_tc = _terminal_vars(m, options)
        weights = _terminal_weights(options, vars_tc, stage_cost_index)

        # Re-detect state_vars after equations_write (user may add DerivativeVar
        # inside equations_write; in practice they never do, but be safe).
        _, state_vars_tc = _get_derivative_and_state_vars(m)

        def _is_differential(var_name):
            base = var_name.split('[')[0]
            return getattr(m, base) in state_vars_tc

        diff_vars = [v for v in vars_tc if _is_differential(v)]
        alg_vars  = [v for v in vars_tc if not _is_differential(v)]

        # Lagrange coefficients: computed once, reused for all alg. vars.
        if alg_vars:
            last_fe_pts = _get_last_fe_pts(m)
            lag_coeffs  = _lagrange_coeffs_at_endpoint(last_fe_pts, 1.0)

            def _alg_extrap(var_name):
                """Pyomo expression for the extrapolated terminal value."""
                return sum(
                    c * _add_time_indexed_expression(m, var_name, t)
                    for c, t in zip(lag_coeffs, last_fe_pts)
                )

        if options.terminal_constraint_type == 'hard':
            # Differential CVs: τ=1 exists and is pinned by Pyomo continuity.
            if diff_vars:
                diff_weights = _terminal_weights(options, diff_vars, stage_cost_index)

                def _diff_tc_rule(m, v):
                    return (
                        _add_time_indexed_expression(m, v, 1)
                        == m.steady_state_values[v]
                    )

                m.diff_terminal_constraints = pyo.Constraint(
                    diff_vars, rule=_diff_tc_rule
                )

            # Algebraic CVs / MVs: polynomial extrapolation to τ=1.
            if alg_vars:
                def _alg_tc_rule(m, v):
                    return _alg_extrap(v) == m.steady_state_values[v]

                m.alg_terminal_constraints = pyo.Constraint(
                    alg_vars, rule=_alg_tc_rule
                )

        else:  # 'soft' — build penalty Expression on the block
            diff_weights = _terminal_weights(options, diff_vars, stage_cost_index)
            alg_weights  = _terminal_weights(options, alg_vars,  stage_cost_index)
            penalty_terms = []
            for v, w in zip(diff_vars, diff_weights):
                penalty_terms.append(
                    w * (_add_time_indexed_expression(m, v, 1)
                         - m.steady_state_values[v])**2
                )
            for v, w in zip(alg_vars, alg_weights):
                penalty_terms.append(
                    w * (_alg_extrap(v) - m.steady_state_values[v])**2
                )
            m.infinite_terminal_soft_penalty = pyo.Expression(
                expr=sum(penalty_terms) if penalty_terms else pyo.Constant(0)
            )

    _transform_model_derivatives(m, deriv_vars)

    return m


def _link_blocks(m):
    """
    Add equality constraints linking state variables at the finite/infinite block interface.

    The last time point of ``finite_block`` is constrained to equal the first
    time point of ``infinite_block`` for each differential state variable.

    Parameters
    ----------
    m : pyo.ConcreteModel

    Returns
    -------
    pyo.ConcreteModel
    """
    t_finite_end = max(m.finite_block.time)
    t_infinite_start = min(m.infinite_block.time)

    def make_link_constraint(name):
        var_finite = getattr(m.finite_block, name)
        var_infinite = getattr(m.infinite_block, name)

        unique_spatial_indices = set()
        for idx in var_finite.index_set():
            if isinstance(idx, tuple):
                unique_spatial_indices.add(idx[:-1])
            else:
                unique_spatial_indices.add(())

        def rule(m, *spatial_idx):
            if spatial_idx:
                return var_finite[spatial_idx + (t_finite_end,)] == \
                       var_infinite[spatial_idx + (t_infinite_start,)]
            else:
                return var_finite[t_finite_end] == var_infinite[t_infinite_start]

        return pyo.Constraint(list(unique_spatial_indices), rule=rule)

    state_var_names = {var.name.split(".")[-1] for var in m.finite_block.state_vars}
    for name in state_var_names:
        setattr(m, f"link_{name}_constraint", make_link_constraint(name))

    return m


# ---- Steady-state economic check ----
if __name__ == '__main__':
    import sys
    from pathlib import Path

    # Make sure the repo root is on sys.path when invoked directly.
    _repo_root = str(Path(__file__).parent.parent)
    if _repo_root not in sys.path:
        sys.path.insert(0, _repo_root)

    from infeNMPC.infNMPC_options import Options  # type: ignore[import]

    # Load the ternary distillation model from the lyap_flag example.
    _example_dir = (
        Path(__file__).parent.parent
        / 'examples' / 'lyap_flag_examples' / 'enmpc_ternary_distillation'
    )
    sys.path.insert(0, str(_example_dir))
    import ternary_distillation_model  # type: ignore[import]  # noqa: E402

    options = Options.for_model_module(
        ternary_distillation_model,
        sampling_time=1,
        nfe_finite=3,
        ncp_finite=3,
        infinite_horizon=False,
        objective='economic',
        tee_flag=True,
        stage_cost_weights=[1.0e4, 1.0e4, 1.0e4, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    )

    print('\nBuilding steady-state model for ternary distillation...')
    m_ss = pyo.ConcreteModel()
    m_ss = _make_steady_state_model(m_ss, options)
    # m_ss.xC.fix(0.95)
    # m_ss.xD1A.fix(0.95)
    # m_ss.xD2B.fix(0.95)

    print('Solving steady-state model...')
    ss_data = _solve_steady_state_model(m_ss, None, options, label="ss")

    # ---- Report objective sign ----
    obj_val = m_ss.ss_obj_value
    print(f'\nSteady-state objective value: {obj_val:.6g}')
    if obj_val < 0:
        print('  -> NEGATIVE: revenue > cost (system is profitable at SS optimum)')
    else:
        print('  -> POSITIVE: cost > revenue (system loses money at SS optimum)')

    # ---- Economic breakdown ----
    # RADAU nfe=1 ncp=1: single collocation point is at t=1.
    t_ss = m_ss.time.last()
    pF, pV, pA, pB, pC = 1.0, 0.008, 1.0, 2.0, 1.0

    F_val   = pyo.value(m_ss.F)
    VB1_val = pyo.value(m_ss.VB1[t_ss])
    VB2_val = pyo.value(m_ss.VB2[t_ss])
    D1_val  = pyo.value(m_ss.D1[t_ss])
    D2_val  = pyo.value(m_ss.D2[t_ss])
    B2_val  = pyo.value(m_ss.B2[t_ss])
    xD1A    = pyo.value(m_ss.xD1A[t_ss])
    xD2B    = pyo.value(m_ss.xD2B[t_ss])
    xC      = pyo.value(m_ss.xC[t_ss])

    feed_cost   = pF * F_val
    energy_cost = pV * (VB1_val + VB2_val)
    revenue     = pA * D1_val + pB * D2_val + pC * B2_val

    print('\nEconomic breakdown at SS optimum (t = t_last):')
    print(f'  Feed cost    pF*F             = {feed_cost:.4f}')
    print(f'  Energy cost  pV*(VB1+VB2)     = {pV}*({VB1_val:.4f} + {VB2_val:.4f}) = {energy_cost:.4f}')
    print(f'  Revenue      pA*D1+pB*D2+pC*B2 = {pA*D1_val:.4f} + {pB*D2_val:.4f} + {pC*B2_val:.4f} = {revenue:.4f}')
    print(f'  Net profit rate (revenue - cost) = {revenue - feed_cost - energy_cost:.4f} $/h')
    print(f'  Objective (cost - revenue)       = {feed_cost + energy_cost - revenue:.4f}')
    print(f'\nPurity CVs: xD1A={xD1A:.4f}  xD2B={xD2B:.4f}  xC={xC:.4f}')
