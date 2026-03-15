"""
Pyomo model construction for infeNMPC.

Provides functions to build steady-state, finite-horizon, and
infinite-horizon NMPC models.  Model equations and variable declarations
are loaded via ``_get_model(options)``, which resolves ``options.model_module``
directly (primary) or falls back to ``_load_model`` by name (legacy).

Non-collocation cleanup is handled by passing ``clean_model='delete'`` to
``dae.collocation.apply_to`` (requires the pyomo-dae ``dae-rewrite`` branch).
``DynamicModelInterface`` is created with ``clean_model=True`` so that all
MPC data-access operations use the safe sparse-access paths from
``contrib.mpc.interfaces.*_clean``.
"""
from .infNMPC_options import _import_settings
import pyomo.environ as pyo
from .model_equations import _get_model
from pyomo.contrib.mpc import DynamicModelInterface, ScalarData
from pyomo.dae import ContinuousSet, DerivativeVar
from pyomo.common.collections import ComponentMap
from pyomo.core.expr.visitor import replace_expressions
from pyomo.core.expr.visitor import identify_variables
from .indexing_tools import (
    _get_variable_key_for_data,
    _add_time_indexed_expression,
    _get_derivative_and_state_vars,
)


def _ipopt_solver():
    """Return a pre-configured IPOPT solver instance."""
    solver = pyo.SolverFactory('ipopt')
    solver.options['linear_solver'] = 'ma57'
    solver.options['OF_ma57_automatic_scaling'] = 'yes'
    solver.options['max_iter'] = 6000
    solver.options['halt_on_ampl_error'] = 'yes'
    solver.options['bound_relax_factor'] = 0
    return solver


def _make_steady_state_model(m, options):
    """
    Generate a steady-state model by fixing all derivatives to zero.

    Discretizes with a single finite element / collocation point so that the
    model represents a static NLP, then deactivates the discretization
    equations so the solver treats derivatives as parameters pinned to zero.

    Parameters
    ----------
    m : pyo.ConcreteModel
    options : Options

    Returns
    -------
    pyo.ConcreteModel
    """
    _model = _get_model(options)

    m.time = ContinuousSet(bounds=(0, 1))
    m = _model.variables_initialize(m)

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


def _solve_steady_state_model(m, target, options):
    """
    Solve the steady-state model and return data at the initial time point.

    When ``options.custom_objective`` is True and no target is supplied, the
    objective is the average of the user-supplied stage cost at ``t=0`` and
    ``t=1``.  Otherwise a standard penalty-from-target objective is used.

    Parameters
    ----------
    m : pyo.ConcreteModel
    target : ScalarData or None
    options : Options

    Returns
    -------
    ScalarData
        Variable values at ``m.time.first()``.
    """
    ss_interface = DynamicModelInterface(m, m.time)

    if options.custom_objective and target is None:
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

    solver = _ipopt_solver()
    solver.solve(m, tee=options.tee_flag)

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

    m_ss_target = None if options.custom_objective else ScalarData({
        _get_variable_key_for_data(m_ss, var): pyo.value(m_ss.setpoints[var])
        for var in m_ss.setpoints.index_set()
    })

    print('Solving Steady State Model')
    steady_state_data = _solve_steady_state_model(m_ss, m_ss_target, options)

    steady_state_values = {
        **{cv: steady_state_data.get_data_from_key(_get_variable_key_for_data(m_ss, cv))
           for cv in m_ss.CV_index},
        **{mv: steady_state_data.get_data_from_key(_get_variable_key_for_data(m_ss, mv))
           for mv in m_ss.MV_index},
    }

    # --- Initial condition solve ---
    print('Writing Initial Value Model')
    m_iv = pyo.ConcreteModel()
    m_iv = _make_steady_state_model(m_iv, options)

    initial_value_vars = list(m_iv.initial_values.index_set())
    if all(var in m_iv.CV_index for var in initial_value_vars):
        initial_data = ScalarData({
            _get_variable_key_for_data(m_iv, cv): pyo.value(m_iv.initial_values[cv])
            for cv in initial_value_vars
        })
    else:
        m_iv_target = ScalarData({
            _get_variable_key_for_data(m_iv, var): pyo.value(m_iv.initial_values[var])
            for var in initial_value_vars
        })
        print('Solving Initial Value Model')
        initial_data = _solve_steady_state_model(m_iv, m_iv_target, options)

    print('Writing Infinite Horizon Model')

    m.finite_block = pyo.Block()
    m.infinite_block = pyo.Block()
    m.finite_block.steady_state_values = steady_state_values
    m.infinite_block.steady_state_values = steady_state_values

    if options.custom_objective:
        m.finite_block.ss_obj_value = m_ss.ss_obj_value
        m.infinite_block.ss_obj_value = m_ss.ss_obj_value

    m.finite_block = _finite_block_gen(m.finite_block, options)
    m.infinite_block = _infinite_block_gen(m.infinite_block, options)
    m = _link_blocks(m)

    m.interface = DynamicModelInterface(
        m.finite_block, m.finite_block.time, clean_model=True
    )

    if options.initialize_with_initial_data:
        for time in m.finite_block.time:
            m.interface.load_data(initial_data, time_points=time)
    else:
        m.interface.load_data(initial_data, time_points=0)

    time0 = m.finite_block.time.first()
    for var in m.finite_block.state_vars:
        for index in var:
            if (isinstance(index, tuple) and index[-1] == time0) or index == time0:
                var[index].fix()

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

    m_ss_target = None if options.custom_objective else ScalarData({
        _get_variable_key_for_data(m_ss, var): pyo.value(m_ss.setpoints[var])
        for var in m_ss.setpoints.index_set()
    })

    steady_state_data = _solve_steady_state_model(m_ss, m_ss_target, options)

    steady_state_values = {
        **{cv: steady_state_data.get_data_from_key(_get_variable_key_for_data(m_ss, cv))
           for cv in m_ss.CV_index},
        **{mv: steady_state_data.get_data_from_key(_get_variable_key_for_data(m_ss, mv))
           for mv in m_ss.MV_index},
    }

    m_iv = pyo.ConcreteModel()
    m_iv = _make_steady_state_model(m_iv, options)

    initial_value_vars = list(m_iv.initial_values.index_set())
    if all(var in m_iv.CV_index for var in initial_value_vars):
        initial_data = ScalarData({
            _get_variable_key_for_data(m_iv, cv): pyo.value(m_iv.initial_values[cv])
            for cv in initial_value_vars
        })
    else:
        m_iv_target = ScalarData({
            _get_variable_key_for_data(m_iv, var): pyo.value(m_iv.initial_values[var])
            for var in initial_value_vars
        })
        initial_data = _solve_steady_state_model(m_iv, m_iv_target, options)

    print('Writing Finite Horizon Model')

    m.ss_obj_value = m_ss.ss_obj_value
    m.steady_state_values = steady_state_values
    m = _finite_block_gen(m, options)

    m.interface = DynamicModelInterface(m, m.time, clean_model=True)

    if options.initialize_with_initial_data:
        for time in m.time:
            m.interface.load_data(initial_data, time_points=time)
    else:
        m.interface.load_data(initial_data, time_points=0)

    time0 = m.time.first()
    for var in m.state_vars:
        for index in var:
            if (isinstance(index, tuple) and index[-1] == time0) or index == time0:
                var[index].fix()

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

    _model.equations_write(m)

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
        if con.name.endswith("_disc_eq") or con.name.endswith("terminal_cost"):
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

    # If no endpoint constraint, remove all entries at t=1
    if not model.endpoint_constraints:
        for comp in model.component_objects(active=True, descend_into=True):
            if not comp.is_indexed():
                continue
            for index in list(comp.keys()):
                time_index = index[-1] if isinstance(index, tuple) else index
                if time_index == 1:
                    try:
                        del comp[index]
                    except (KeyError, AttributeError, TypeError):
                        pass


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

    def _add_dphidt(m):
        m.phi = pyo.Var(m.time, initialize=0, domain=pyo.NonNegativeReals)
        m.phi[0].fix(0)
        m.dphidt = DerivativeVar(m.phi, wrt=m.time)

        m.stage_cost_index = pyo.Set(
            initialize=lambda m: list(m.CV_index) + list(m.MV_index)
        )
        c = options.stage_cost_weights

        def _terminal_cost_rule(m, t):
            stage_cost = sum(
                c[i] * (
                    _add_time_indexed_expression(m, var_name, t)
                    - m.steady_state_values[var_name]
                )**2
                for i, var_name in enumerate(m.stage_cost_index)
            )
            if t == 1 and options.endpoint_constraints:
                if any(c[i] == 0 for i in range(len(c))):
                    cmod = [1 for _ in c]
                    stage_cost_mod = sum(
                        cmod[i] * (
                            _add_time_indexed_expression(m, var_name, t)
                            - m.steady_state_values[var_name]
                        )**2
                        for i, var_name in enumerate(m.stage_cost_index)
                    )
                    return stage_cost_mod <= 1e-6
                else:
                    return stage_cost <= 1e-6
            elif t == 0:
                return pyo.Constraint.Skip
            else:
                return (
                    (options.gamma / options.sampling_time * (1 - t**2))
                    * m.dphidt[t] == stage_cost
                )

        def _custom_terminal_cost_rule(m, t):
            setpoint_stage_cost = sum(
                c[i] * (
                    _add_time_indexed_expression(m, var_name, t)
                    - m.steady_state_values[var_name]
                )**2
                for i, var_name in enumerate(m.stage_cost_index)
            )
            # custom_objective is bound via closure from the if-block below
            cost_fn = custom_objective(m, options)
            if t == 1:
                if options.endpoint_constraints:
                    return setpoint_stage_cost <= 1e-6
                else:
                    return pyo.Constraint.Skip
            elif t == 0:
                return pyo.Constraint.Skip
            else:
                stage_cost_expr = cost_fn(m, t)
                return (
                    (options.gamma / options.sampling_time * (1 - t**2))
                    * m.dphidt[t] == stage_cost_expr - m.ss_obj_value
                )

        if options.custom_objective:
            custom_objective = _model.custom_objective
            m.terminal_cost = pyo.Constraint(
                m.time, rule=_custom_terminal_cost_rule
            )
        else:
            m.terminal_cost = pyo.Constraint(m.time, rule=_terminal_cost_rule)

        return m

    deriv_vars, state_vars = _get_derivative_and_state_vars(m)
    m = _add_dphidt(m)

    m.gamma = options.gamma
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

    m.endpoint_constraints = options.endpoint_constraints
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


# ---- Testing the Model ----
if __name__ == '__main__':
    options = _import_settings()
    m = pyo.ConcreteModel()

    print('\nTesting Finite Horizon Model Setup')
    m = _make_finite_horizon_model(m, options)

    if True:
        with open("model_output.txt", "w") as f:
            m.pprint(ostream=f)

    assert isinstance(m, pyo.ConcreteModel)
