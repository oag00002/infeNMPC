"""
Progressive sampling-time warm-start utilities for robust controller initialisation.
"""
from ..infNMPC_options import Options
from pyomo.contrib.mpc.data.series_data import TimeSeriesData


def _assist_initialization_infinite(options: Options):
    """
    Warm-start an infinite-horizon controller via progressive sampling-time reduction.

    Solves a sequence of controllers beginning at
    ``options.initialization_assist_sampling_time_start`` and iteratively
    multiplying the sampling time by 0.9 until the target
    ``options.sampling_time`` is reached.  Each solution is loaded as the
    initial guess for the next, improving convergence for stiff problems.

    Parameters
    ----------
    options : Options
        Must include ``initialization_assist_sampling_time_start`` in addition
        to all standard controller options.

    Returns
    -------
    InfiniteHorizonController
        The solved controller at the target sampling time.
    """
    from ..controllers import InfiniteHorizonController

    print("Using initialization assist for infinite horizon controller.")

    modified_options = options.copy(
        sampling_time=options.initialization_assist_sampling_time_start
    )

    # Create initial controller and extract data
    controller = InfiniteHorizonController(modified_options)
    interface = controller.interface
    time_points = (
        list(controller.finite_block.time) + list(controller.infinite_block.time)
    )
    print("Initial time points:", time_points)
    previous_time = list(controller.finite_block.time) + list(controller.infinite_block.time)
    solution_data_list = [interface.get_data_at_time(t) for t in previous_time]

    # Iteratively reduce sampling time
    while modified_options.sampling_time * 0.9 > options.sampling_time:
        modified_options = modified_options.copy(
            sampling_time=round(modified_options.sampling_time * 0.9, 8)
        )
        print(f"  Trying sampling time: {modified_options.sampling_time}")

        controller = InfiniteHorizonController(modified_options)
        interface = controller.interface
        current_time = (
            list(controller.finite_block.time) + list(controller.infinite_block.time)
        )
        print("Updated time points:", current_time)
        new_data = TimeSeriesData(
            {cuid: [d[cuid] for d in solution_data_list] for cuid in solution_data_list[0]},
            time=current_time,
            time_set=interface.time,
        )
        interface.load_data(new_data, time_points=current_time)
        solution_data_list = [
            interface.get_data_at_time(t)
            for t in list(controller.finite_block.time) + list(controller.infinite_block.time)
        ]

    # Final solve at target sampling time
    print(f"Final solve at target sampling time: {options.sampling_time}")
    controller = InfiniteHorizonController(options)
    interface = controller.interface
    final_time = (
        list(controller.finite_block.time) + list(controller.infinite_block.time)
    )
    final_data = TimeSeriesData(
        {cuid: [d[cuid] for d in solution_data_list] for cuid in solution_data_list[0]},
        time=final_time,
        time_set=interface.time,
    )
    interface.load_data(final_data, time_points=final_time)
    return controller


def _assist_initialization_finite(options: Options):
    """
    Warm-start a finite-horizon controller via progressive sampling-time reduction.

    Mirrors ``_assist_initialization_infinite`` for the finite-horizon case.

    Parameters
    ----------
    options : Options
        Must include ``initialization_assist_sampling_time_start`` in addition
        to all standard controller options.

    Returns
    -------
    FiniteHorizonController
        The solved controller at the target sampling time.
    """
    from ..controllers import FiniteHorizonController

    print("Using initialization assist for finite horizon controller.")

    modified_options = options.copy(
        sampling_time=options.initialization_assist_sampling_time_start
    )

    # Create initial controller and extract data
    controller = FiniteHorizonController(modified_options)
    previous_time = list(controller.time)
    print("Initial time points:", previous_time)
    solution_data_list = [controller.interface.get_data_at_time(t) for t in previous_time]

    # Iteratively reduce sampling time
    while modified_options.sampling_time * 0.9 > options.sampling_time:
        modified_options = modified_options.copy(
            sampling_time=round(modified_options.sampling_time * 0.9, 8)
        )
        print(f"  Trying sampling time: {modified_options.sampling_time}")

        controller = FiniteHorizonController(modified_options, solution_data_list)
        current_time = list(controller.time)
        print("Updated time points:", current_time)
        solution_data_list = [controller.interface.get_data_at_time(t) for t in current_time]

    # Final solve at target sampling time
    print(f"Final solve at target sampling time: {options.sampling_time}")
    return FiniteHorizonController(options, solution_data_list)


def _check_ic_consistency(block, tol=1e-4):
    """
    Evaluate all active equality constraints at t=0 after ICs are fixed.

    Must be called after ``_fix_initial_conditions_at_t0`` has fixed all
    differential state variables at t=0.  Algebraic variables at t=0 have
    been loaded (not fixed) from ``initial_data``; ``pyo.value()`` reads
    their current numerical values.

    Discretization constraints (names ending in ``_disc_eq``) are skipped
    because they encode the collocation rule, not the physical equations.

    Parameters
    ----------
    block : pyo.ConcreteModel or pyo.Block
        A discretised model block with ``.time``.
    tol : float
        Residual magnitude above which a constraint is considered violated.

    Raises
    ------
    RuntimeError
        If any active equality constraint at t=0 has |residual| > tol.
        The message lists the top 5 violations by magnitude.
    """
    import pyomo.environ as pyo

    t0 = block.time.first()
    violations = []

    for con in block.component_objects(pyo.Constraint, active=True):
        if con.name.endswith("_disc_eq"):
            continue
        for idx in con:
            time_val = idx[-1] if isinstance(idx, tuple) else idx
            if time_val != t0:
                continue
            con_data = con[idx]
            lower = pyo.value(con_data.lower) if con_data.lower is not None else None
            upper = pyo.value(con_data.upper) if con_data.upper is not None else None
            if lower is None or upper is None or lower != upper:
                continue  # skip inequalities
            try:
                body_val = pyo.value(con_data.body, exception=False)
            except Exception:
                continue
            if body_val is None:
                continue
            residual = abs(body_val - lower)
            if residual > tol:
                violations.append((residual, con.name, idx))

    if violations:
        violations.sort(key=lambda x: -x[0])
        lines = [f"  {res:.6g}  {name}[{idx}]" for res, name, idx in violations[:5]]
        raise RuntimeError(
            "Initial conditions are inconsistent at t=0. "
            "Top constraint violations:\n" + "\n".join(lines)
        )


def _build_ic_data_steady_state(options):
    """
    Build initial condition data by solving a steady-state NLP toward
    ``initial_values``, with slack variables relaxed (``fix_slacks=False``).

    This is the default IC path: IPOPT finds a consistent static point near
    the user-supplied ``initial_values``.  Because slacks are free, the IV
    model can converge even when the target initial point violates soft
    constraints (e.g., product purity specs in distillation).

    Parameters
    ----------
    options : Options

    Returns
    -------
    ScalarData
        Variable values at t_first after a successful solve.

    Raises
    ------
    RuntimeError
        With a user-readable message if IPOPT fails.
    """
    import pyomo.environ as pyo
    from pyomo.contrib.mpc import ScalarData
    from ..make_model import _make_steady_state_model, _solve_steady_state_model
    from .indexing_tools import _get_variable_key_for_data

    print('Writing Initial Value Model')
    m_iv = pyo.ConcreteModel()
    m_iv = _make_steady_state_model(m_iv, options, fix_slacks=False)

    initial_value_vars = list(m_iv.initial_values.index_set())
    m_iv_target = ScalarData({
        _get_variable_key_for_data(m_iv, var): pyo.value(m_iv.initial_values[var])
        for var in initial_value_vars
    })

    print('Solving Initial Value Model')
    try:
        initial_data = _solve_steady_state_model(m_iv, m_iv_target, options, label="iv")
    except RuntimeError as e:
        raise RuntimeError(
            "Initial condition solve failed — the specified initial_values may define "
            "an infeasible or inconsistent starting point. Check that initial_values "
            "are achievable under the model equations.\n"
            f"Original: {e}"
        ) from e

    return initial_data


def _build_ic_data_dynamic(options):
    """
    Build initial condition data for a non-steady-state starting point.

    State variables are fixed to their ``initialize=`` values (from
    ``variables_initialize``), with optional scalar overrides from
    ``initial_values``.  Algebraic variables and derivatives are then solved
    for via a feasibility problem: ``disc_eq`` is deactivated so derivatives
    are free variables determined by the ODE constraints at t=0.  This allows
    the starting point to be off steady-state.

    Works for both scalar and spatially-indexed state variables (e.g.,
    distillation tray compositions) — for indexed vars, the ``initialize=``
    values set in ``variables_initialize`` encode the per-index IC.
    ``initial_values`` only provides scalar overrides for specific variables
    (e.g., CV targets such as ``xD1A=0.90``).

    Parameters
    ----------
    options : Options
        Must have ``options.dynamic_initial_conditions = True``.

    Returns
    -------
    ScalarData
        All variable values at t_first after a successful feasibility solve.

    Raises
    ------
    RuntimeError
        With constraint violation details if IPOPT cannot satisfy the
        algebraic constraints at the specified initial state.
    """
    import pyomo.environ as pyo
    from pyomo.dae import ContinuousSet
    from pyomo.contrib.mpc import DynamicModelInterface
    from pyomo.opt import TerminationCondition
    from ..model_equations import _get_model
    from .indexing_tools import _get_derivative_and_state_vars
    from .debug_tools import _report_constraint_violations
    from ..make_model import _ipopt_solver

    print('Writing Dynamic IC Model (non-steady-state specification)')
    _model = _get_model(options)

    m_iv = pyo.ConcreteModel()
    m_iv.time = ContinuousSet(bounds=(0, 1))
    m_iv = _model.variables_initialize(m_iv)

    deriv_vars, state_var_components = _get_derivative_and_state_vars(m_iv)

    # Apply any scalar overrides from initial_values (e.g., CV targets)
    for iv_name in m_iv.initial_values.index_set():
        iv_val = pyo.value(m_iv.initial_values[iv_name])
        var = getattr(m_iv, iv_name, None)
        if var is not None and isinstance(var, pyo.Var):
            for idx in var.index_set():
                var[idx].set_value(iv_val)

    # Discretize (nfe=1, ncp=1) — same mesh as steady-state model
    discretizer = pyo.TransformationFactory('dae.collocation')
    discretizer.apply_to(m_iv, ncp=1, nfe=1, wrt=m_iv.time, scheme='LAGRANGE-RADAU')

    _model.equations_write(m_iv)

    # Deactivate disc_eq — derivatives become free, determined by ODE constraints
    for var in deriv_vars:
        disc_eq_name = f"{var.name.split('.')[-1]}_disc_eq"
        disc_eq = getattr(m_iv, disc_eq_name, None)
        if disc_eq is not None:
            disc_eq.deactivate()

    t_first = m_iv.time.first()

    # Fix all differential state vars at t_first to current (initialize=) values
    for sv in state_var_components:
        for idx in sv.index_set():
            time_val = idx[-1] if isinstance(idx, tuple) else idx
            if time_val == t_first:
                current_val = pyo.value(sv[idx])
                if current_val is not None:
                    sv[idx].fix(current_val)

    # Fix MVs to prevent underdetermined system
    for mv_name in m_iv.MV_index:
        mv_var = getattr(m_iv, mv_name)
        for idx in mv_var.index_set():
            time_val = idx[-1] if isinstance(idx, tuple) else idx
            if time_val == t_first:
                current_val = pyo.value(mv_var[idx])
                if current_val is None and mv_var[idx].lb is not None:
                    current_val = mv_var[idx].lb
                mv_var[idx].fix(current_val if current_val is not None else 0.0)

    # Feasibility objective — IPOPT solves for algebraic vars and free derivatives
    m_iv.objective = pyo.Objective(expr=0)

    if options.model_output_dir:
        import os
        os.makedirs(options.model_output_dir, exist_ok=True)
        with open(os.path.join(options.model_output_dir, "dynamic_ic_model.txt"), "w") as f:
            m_iv.pprint(ostream=f)

    solver = _ipopt_solver()
    results = solver.solve(m_iv, tee=options.tee_flag)
    tc = results.solver.termination_condition

    if tc != TerminationCondition.optimal:
        _report_constraint_violations(m_iv, label="dynamic IC infeasibility")
        raise RuntimeError(
            "Dynamic initial conditions are inconsistent: the specified state "
            "variable values violate algebraic constraints. "
            "Check that your initialize= values are physically realizable. "
            f"Solver termination: {tc}"
        )

    iv_interface = DynamicModelInterface(m_iv, m_iv.time)
    return iv_interface.get_data_at_time(t_first)


def _build_ic_data(options):
    """
    Route to the appropriate initial condition data builder.

    Calls ``_build_ic_data_dynamic`` when
    ``options.dynamic_initial_conditions`` is True, otherwise calls
    ``_build_ic_data_steady_state``.

    Parameters
    ----------
    options : Options

    Returns
    -------
    ScalarData
        Variable values to load at t=0 on the finite block.
    """
    if options.dynamic_initial_conditions:
        return _build_ic_data_dynamic(options)
    else:
        return _build_ic_data_steady_state(options)
