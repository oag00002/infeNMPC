from infNMPC_options import _import_settings
import pyomo.environ as pyo
from model_equations import equations_write, variables_initialize
import idaes
from idaes.core.solvers import use_idaes_solver_configuration_defaults
from pyomo.contrib.mpc import DynamicModelInterface, ScalarData
from pyomo.dae import ContinuousSet, DerivativeVar
from pyomo.common.collections import ComponentMap
from pyomo.core.expr.visitor import replace_expressions
from pyomo.core.expr.visitor import identify_variables
from indexing_tools import _get_variable_key_for_data, _add_time_indexed_expression, _get_derivative_and_state_vars, _get_disc_eq_time_points


# Solver Settings
use_idaes_solver_configuration_defaults()
idaes.cfg.ipopt.options.linear_solver = "ma57" # "ma27"
idaes.cfg.ipopt.options.OF_ma57_automatic_scaling = "yes"
idaes.cfg.ipopt.options.max_iter = 6000
idaes.cfg.ipopt.options.halt_on_ampl_error = "yes"
idaes.cfg.ipopt.options.bound_relax_factor = 0


def _make_steady_state_model(m, options):
    """ 
    Generate steady state values for model
    
    Parameters
    ----------
    m : ConcreteModel
        The Pyomo model object to be modified.
        
    Returns
    -------
    ConcreteModel
        The modified Pyomo model with steady state values.
    """
    m.time = ContinuousSet(bounds=(0, 1))

    m = variables_initialize(m)

    deriv_vars, state_vars = _get_derivative_and_state_vars(m)

    discretizer = pyo.TransformationFactory('dae.collocation')
    discretizer.apply_to(m, ncp=1, nfe=1, wrt=m.time, scheme='LAGRANGE-RADAU')

    m = equations_write(m)

    for var in deriv_vars:
        for idx in var.index_set():
            # Check if last index is time
            if isinstance(idx, tuple):
                time_val = idx[-1]
            else:
                time_val = idx

            if time_val in m.time:
                var[idx].fix(0)

    for var in deriv_vars:
        getattr(m, f"{var}_disc_eq").deactivate()

    return m


def _solve_steady_state_model(m, target, options):
    """
    Solve the steady state model to obtain steady state values.

    Parameters
    ----------
    m : ConcreteModel
        The Pyomo model object to be solved.

    Returns
    -------
    ConcreteModel
        The modified Pyomo model with solved steady state values.
    """
    ss_interface = DynamicModelInterface(m, m.time)

    if options.custom_objective and target is None:
        from model_equations import custom_objective
        cost_fn = custom_objective(m, options)

        def ss_obj_rule(m):
            return (cost_fn(m, 0) + cost_fn(m, 1)) / 2
        
        m.objective = pyo.Objective(expr=ss_obj_rule(m), sense=pyo.minimize)

    else:
        var_set, tr_cost = ss_interface.get_penalty_from_target(target)

        m.target_set = var_set
        m.tracking_cost = tr_cost
        m.objective = pyo.Objective(expr=sum(m.tracking_cost[:, 0]))
    
    solver = pyo.SolverFactory('ipopt')
    solver.solve(m, tee=options.tee_flag)

    m.ss_obj_value = pyo.value(m.objective)

    steady_state_data = ss_interface.get_data_at_time(m.time.first())

    # Extract MV and CV values from solved steady state
    # steady_state_input_data = {
    #     mv: steady_state_data.get_data_from_key(_get_variable_key_for_data(m, mv))
    #     for mv in m.MV_index
    # }
    # steady_state_output_data = {
    #     cv: steady_state_data.get_data_from_key(_get_variable_key_for_data(m, cv))
    #     for cv in m.CV_index
    # }
    # steady_state_values = {**steady_state_output_data, **steady_state_input_data}

    # print(steady_state_values)

    return steady_state_data


def _make_infinite_horizon_model(m, options):
    """
    Create an infinite horizon NMPC model.

    Parameters
    ----------
    m : ConcreteModel
        The Pyomo model object to be modified.

    Returns
    -------
    ConcreteModel
        The modified Pyomo model with infinite horizon NMPC settings.
    """
    # --- Build and solve steady-state model with setpoints ---
    print('Writing Steady State Model')
    m_ss = pyo.ConcreteModel()
    m_ss = _make_steady_state_model(m_ss, options)

    if options.custom_objective:
        m_ss_target = None
    else:
        setpoint_targets = {
            _get_variable_key_for_data(m_ss, var): pyo.value(m_ss.setpoints[var])
            for var in m_ss.setpoints.index_set()
        }
        m_ss_target = ScalarData(setpoint_targets)

    print('Solving Steady State Model')
    steady_state_data = _solve_steady_state_model(m_ss, m_ss_target, options)

    # Extract MV and CV values from solved steady state
    steady_state_input_data = {
        mv: steady_state_data.get_data_from_key(_get_variable_key_for_data(m_ss, mv))
        for mv in m_ss.MV_index
    }
    steady_state_output_data = {
        cv: steady_state_data.get_data_from_key(_get_variable_key_for_data(m_ss, cv))
        for cv in m_ss.CV_index
    }
    steady_state_values = {**steady_state_output_data, **steady_state_input_data}

    print('Writing Initial Value Model')
    # --- Handle initial values ---
    m_iv = pyo.ConcreteModel()
    m_iv = _make_steady_state_model(m_iv, options)

    initial_value_vars = list(m_iv.initial_values.index_set())
    if all(var in m_iv.CV_index for var in initial_value_vars):
        # If CVs are passed as initial values, use CV initial values only
        initial_data_dict = {
            _get_variable_key_for_data(m_iv, cv): pyo.value(m_iv.initial_values[cv])
            for cv in initial_value_vars
        }
        initial_data = ScalarData(initial_data_dict)

    else:
        # Otherwise assume MVs, and solve for initial state
        iv_targets = {
            _get_variable_key_for_data(m_iv, var): pyo.value(m_iv.initial_values[var])
            for var in initial_value_vars
        }
        m_iv_target = ScalarData(iv_targets)
        print('Solving Initial Value Model')
        initial_data = _solve_steady_state_model(m_iv, m_iv_target, options)

    print('Writing Infinite Horizon Model')

    # Define blocks
    m.finite_block = pyo.Block()
    m.infinite_block = pyo.Block()
    m.finite_block.steady_state_values = steady_state_values
    m.infinite_block.steady_state_values = steady_state_values

    if options.custom_objective:
        m.finite_block.ss_obj_value = m_ss.ss_obj_value
        m.infinite_block.ss_obj_value = m_ss.ss_obj_value

    # Generate blocks
    m.finite_block = _finite_block_gen(m.finite_block, options)
    m.infinite_block = _infinite_block_gen(m.infinite_block, options)
    m = _link_blocks(m)

    m.interface = DynamicModelInterface(m.finite_block, m.finite_block.time)
    if options.initialize_with_initial_data:
        for time in m.finite_block.time:
            m.interface.load_data(initial_data, time_points=time)
    else:
        m.interface.load_data(initial_data, time_points=0)

    time0 = m.finite_block.time.first()

    if options.remove_collocation:
        m = _remove_non_collocation_values_infinite(m)
    else:
        for mv in m.MV_index:
            getattr(m,mv)[time0].fix()

    for var in m.finite_block.state_vars:
        for index in var:
            if (isinstance(index, tuple) and index[-1] == time0) or index == time0:
                var[index].fix()

    return m


def _remove_non_collocation_values_infinite(m):

    # ---- Finite Block Pruning ----

    # Prune variables and constraints based on collocation time points
    keep_times_finite = set(_get_disc_eq_time_points(m.finite_block))
    state_var_names = {v.parent_component().name.split('.')[-1] for v in m.finite_block.state_vars}
    # print("State variable component names:", state_var_names)

    for comp in m.finite_block.component_objects(pyo.Var, descend_into=True, active=True):
        comp_base_name = comp.name.split('.')[-1]
        # print("Checking", comp.name)
        if comp_base_name not in state_var_names:
            # print(f"  -> {comp.name} marked for pruning.")
            for idx in list(comp):
                time_val = idx[-1] if isinstance(idx, tuple) else idx
                if time_val in m.finite_block.time and time_val not in keep_times_finite:
                    # print(f"    -> Deleting {comp.name}[{idx}] at time {time_val}")
                    del comp[idx]
        # else:
            # print(f"  -> {comp.name} is a state variable component. Skipping.")


    for con in m.finite_block.component_objects(pyo.Constraint, active=True):
        if not con.name.endswith('_time_cont_eq'):
            for idx in list(con):
                time_val = idx[-1] if isinstance(idx, tuple) else idx
                if time_val in m.finite_block.time and time_val not in keep_times_finite:
                    del con[idx]

    # ---- Infinite Block Pruning ----

    # Prune variables and constraints based on collocation time points
    keep_times_infinite = set(_get_disc_eq_time_points(m.infinite_block))
    # print(f"Keeping Times {keep_times_infinite}")
    # print("State variable component names:", state_var_names)

    for comp in m.infinite_block.component_objects(pyo.Var, descend_into=True, active=True):
        comp_base_name = comp.name.split('.')[-1]
        # print("Checking", comp.name)
        if comp_base_name not in state_var_names and comp_base_name != "phi":
            # print(f"  -> {comp.name} marked for pruning.")
            for idx in list(comp):
                # print(f"   -> Index is {idx}")
                time_val = idx[-1] if isinstance(idx, tuple) else idx
                if time_val in m.infinite_block.time and time_val not in keep_times_infinite:
                    # print(f"    -> Deleting {comp.name}[{idx}] at time {time_val}")
                    del comp[idx]
        # else:
            # print(f"  -> {comp.name} is a state variable component. Skipping.")

    for con in m.infinite_block.component_objects(pyo.Constraint, active=True):
        # Skip time continuity constraints
        if con.name.endswith('_time_cont_eq'):
            continue

        for idx in list(con):
            # Extract time value from index
            time_val = idx[-1] if isinstance(idx, tuple) else idx

            # Ignore if it's not a time index
            if time_val not in m.infinite_block.time:
                continue

            if con.name.endswith('terminal_cost'):
                # Keep if time == 1 or time is in keep_times_infinite
                if time_val == 1 or time_val in keep_times_infinite:
                    continue  # DO NOT delete
                else:
                    del con[idx]  # DELETE this terminal cost
            else:
                # For all other constraints, delete if time not in keep_times_infinite
                if time_val not in keep_times_infinite:
                    del con[idx]

    
    return m


def _remove_non_collocation_values_finite(m):

    # ---- Finite Block Pruning ----

    # Prune variables and constraints based on collocation time points
    keep_times_finite = set(_get_disc_eq_time_points(m))
    state_var_names = {v.name.split('.')[-1] for v in m.state_vars}
    # print("State variable component names:", state_var_names)

    for comp in m.component_objects(pyo.Var, descend_into=True, active=True):
        comp_base_name = comp.name.split('.')[-1]
        # print("Checking", comp.name)
        if comp_base_name not in state_var_names:
            # print(f"  -> {comp.name} marked for pruning.")
            for idx in list(comp):
                time_val = idx[-1] if isinstance(idx, tuple) else idx
                if time_val in m.time and time_val not in keep_times_finite:
                    # print(f"    -> Deleting {comp.name}[{idx}] at time {time_val}")
                    del comp[idx]
        # else:
            # print(f"  -> {comp.name} is a state variable component. Skipping.")


    for con in m.component_objects(pyo.Constraint, active=True):
        if not con.name.endswith('_time_cont_eq'):
            for idx in list(con):
                time_val = idx[-1] if isinstance(idx, tuple) else idx
                if time_val in m.time and time_val not in keep_times_finite:
                    del con[idx]
    
    return m


def _make_finite_horizon_model(m, options):
    """
    Create a finite horizon NMPC model.

    Parameters
    ----------
    m : ConcreteModel
        The Pyomo model object to be modified.

    Returns
    -------
    ConcreteModel
        The modified Pyomo model with infinite horizon NMPC settings.
    """
    # Create steady-state model

    # --- Build and solve steady-state model with setpoints ---
    m_ss = pyo.ConcreteModel()
    m_ss = _make_steady_state_model(m_ss, options)

    if options.custom_objective:
        m_ss_target = None
    else:
        setpoint_targets = {
            _get_variable_key_for_data(m_ss, var): pyo.value(m_ss.setpoints[var])
            for var in m_ss.setpoints.index_set()
        }
        # print(setpoint_targets)
        m_ss_target = ScalarData(setpoint_targets)

    steady_state_data = _solve_steady_state_model(m_ss, m_ss_target, options)

    # Extract MV and CV values from solved steady state
    steady_state_input_data = {
        mv: steady_state_data.get_data_from_key(_get_variable_key_for_data(m_ss, mv))
        for mv in m_ss.MV_index
    }
    steady_state_output_data = {
        cv: steady_state_data.get_data_from_key(_get_variable_key_for_data(m_ss, cv))
        for cv in m_ss.CV_index
    }
    steady_state_values = {**steady_state_output_data, **steady_state_input_data}

    # --- Handle initial values ---
    m_iv = pyo.ConcreteModel()
    m_iv = _make_steady_state_model(m_iv, options)

    initial_value_vars = list(m_iv.initial_values.index_set())
    if all(var in m_iv.CV_index for var in initial_value_vars):
        # If CVs are passed as initial values, use CV initial values only
        initial_data_dict = {
            _get_variable_key_for_data(m_iv, cv): pyo.value(m_iv.initial_values[cv])
            for cv in initial_value_vars
        }
        initial_data = ScalarData(initial_data_dict)

    else:
        # Otherwise assume MVs, and solve for initial state
        iv_targets = {
            _get_variable_key_for_data(m_iv, var): pyo.value(m_iv.initial_values[var])
            for var in initial_value_vars
        }
        m_iv_target = ScalarData(iv_targets)
        initial_data = _solve_steady_state_model(m_iv, m_iv_target, options)

    print('Writing Finite Horizon Model')

    m.ss_obj_value = m_ss.ss_obj_value

    m.steady_state_values = steady_state_values

    # Generate blocks
    m = _finite_block_gen(m, options)

    m.interface = DynamicModelInterface(m, m.time)
    
    if options.initialize_with_initial_data:
        for time in m.time:
            m.interface.load_data(initial_data, time_points=time)
    else:
        m.interface.load_data(initial_data, time_points=0)

    time0 = m.time.first()

    if options.remove_collocation:
        m = _remove_non_collocation_values_finite(m)
    else:
        for mv in m.MV_index:
            getattr(m,mv)[time0].fix()

    for var in m.state_vars:
        for index in var:
            if (isinstance(index, tuple) and index[-1] == time0) or index == time0:
                var[index].fix()

    return m


def _finite_block_gen(m, options):
    """
    Generate the finite block for the NMPC model.

    Parameters
    ----------
    m : ConcreteModel
        The Pyomo model object to be modified.

    Returns
    -------
    None
        The function modifies the model in place.
    """
    m.time = ContinuousSet(bounds=(0, options.nfe_finite*options.sampling_time))

    m = variables_initialize(m)

    deriv_vars, state_vars = _get_derivative_and_state_vars(m)
    m.deriv_vars = deriv_vars
    m.state_vars = state_vars

    discretizer = pyo.TransformationFactory('dae.collocation')
    discretizer.apply_to(m, ncp=options.ncp_finite, nfe=options.nfe_finite, wrt=m.time, scheme='LAGRANGE-RADAU')
    if options.ncp_finite > 1:
        for mv in m.MV_index:
            m = discretizer.reduce_collocation_points(m, var=getattr(m, mv), ncp=1, contset=m.time)

    equations_write(m)

    return m


def _transform_model_derivatives(model, deriv_vars):
    """
    Transforms constraints by replacing dxdt[...] with gamma*(1 - t**2)*dxdt[...].
    Uses derivative variable names collected before discretization.
    Skips constraints whose names end with '_disc_eq'.
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

            # Safely extract the time index
            if isinstance(index, tuple):
                if len(index) == 0:
                    continue  # skip empty tuple
                time_index = index[-1]
                index_prefix = index[:-1]  # Use all but last for new constraint index
            else:
                time_index = index
                index_prefix = ()  # Use empty tuple for scalar indices

            # Always make index_prefix a tuple
            if not isinstance(index_prefix, tuple):
                index_prefix = (index_prefix,)

            replacement_map = ComponentMap()
            for var in identify_variables(expr, include_fixed=False):
                if var.parent_component() in deriv_vars:
                    replacement_map[var] = gamma / sampling_time * (1 - time_index**2) * var

            if replacement_map:
                new_expr = replace_expressions(expr, substitution_map=replacement_map)
                con[index].set_value(new_expr)
    
    # Delete all components with time_index == 1 if endpoint_constraints is False
    if not model.endpoint_constraints:
        for comp in model.component_objects(active=True, descend_into=True):
            if not comp.is_indexed():
                continue

            for index in list(comp.keys()):
                if isinstance(index, tuple):
                    if len(index) == 0:
                        continue
                    time_index = index[-1]
                else:
                    time_index = index

                if time_index == 1:
                    try:
                        del comp[index]
                    except (KeyError, AttributeError, TypeError):
                        pass  # Component may not support deletion or index doesn't exist


def _infinite_block_gen(m, options):
    m.time = ContinuousSet(bounds=(0,1))

    m = variables_initialize(m)

    def _add_dphidt(m):

        m.phi = pyo.Var(m.time, initialize=0, domain=pyo.NonNegativeReals)
        m.phi[0].fix(0)  # Initial condition for phi
        m.dphidt = DerivativeVar(m.phi, wrt=m.time)

        m.stage_cost_index = pyo.Set(initialize=lambda m: list(m.CV_index) + list(m.MV_index))

        c = options.stage_cost_weights

        def _terminal_cost_rule(m, t):
            stage_cost = sum(
                c[i] * (_add_time_indexed_expression(m, var_name, t) - m.steady_state_values[var_name])**2
                for i, var_name in enumerate(m.stage_cost_index)
            )

            if t == 1 and options.endpoint_constraints:
                if any(c[i] == 0 for i in range(len(c))):
                    cmod = [1 for _ in c]  # Set all cmod[i] = 1
                    stage_cost_mod = sum(
                        cmod[i] * (_add_time_indexed_expression(m, var_name, t) - m.steady_state_values[var_name])**2
                        for i, var_name in enumerate(m.stage_cost_index)
                    )
                    return stage_cost_mod <= 1e-6
                else:
                    return stage_cost <= 1e-6
                # else:
                #     return stage_cost <= 1e-6
            elif t == 0:
                return pyo.Constraint.Skip
            else:
                return (options.gamma / options.sampling_time * (1 - t**2)) * m.dphidt[t] == stage_cost
            
        def _custom_terminal_cost_rule(m, t):
            setpoint_stage_cost = sum(
                c[i] * (_add_time_indexed_expression(m, var_name, t) - m.steady_state_values[var_name])**2
                for i, var_name in enumerate(m.stage_cost_index)
            )
            cost_fn = custom_objective(m, options)  # this returns a lambda
            if t == 1:
                if options.endpoint_constraints:
                    return setpoint_stage_cost <= 1e-6
                else:
                    return pyo.Constraint.Skip
            elif t == 0:
                return pyo.Constraint.Skip
            else:
                stage_cost_expr = cost_fn(m, t)  # now it's an expression
                return (options.gamma / options.sampling_time * (1 - t**2)) * m.dphidt[t] == stage_cost_expr - m.ss_obj_value

        if options.custom_objective:
            from model_equations import custom_objective
            m.terminal_cost = pyo.Constraint(m.time, rule=_custom_terminal_cost_rule)
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
    discretizer.apply_to(m, ncp=options.ncp_infinite, nfe=options.nfe_infinite, wrt=m.time, scheme='LAGRANGE-LEGENDRE')

    equations_write(m)

    m.endpoint_constraints = options.endpoint_constraints

    _transform_model_derivatives(m, deriv_vars)

    return m



def _link_blocks(m):
    """
    Link the differential variables at the interface between the finite and infinite blocks.

    Parameters
    ----------
    m : ConcreteModel
        The Pyomo model object with finite_block and infinite_block.

    Returns
    -------
    None
        Modifies the model in place by adding interface constraints.
    """
    t_finite_end = max(m.finite_block.time)
    t_infinite_start = min(m.infinite_block.time)

    def make_link_constraint(name):
        var_finite = getattr(m.finite_block, name)
        var_infinite = getattr(m.infinite_block, name)

        # Extract spatial indices (excluding time)
        unique_spatial_indices = set()
        for idx in var_finite.index_set():
            if isinstance(idx, tuple):
                unique_spatial_indices.add(idx[:-1])
            else:
                # Scalar variable over time only
                unique_spatial_indices.add(())

        def rule(m, *spatial_idx):
            if spatial_idx:
                finite_index = spatial_idx + (t_finite_end,)
                infinite_index = spatial_idx + (t_infinite_start,)
            else:
                finite_index = t_finite_end
                infinite_index = t_infinite_start
            return var_finite[finite_index] == var_infinite[infinite_index]

        return pyo.Constraint(list(unique_spatial_indices), rule=rule)
    
    state_var_names = {var.name.split(".")[-1] for var in m.finite_block.state_vars}

    for name in state_var_names:
        setattr(m, f"link_{name}_constraint", make_link_constraint(name))

    return m


# ---- Testing the Model ----
if __name__ == '__main__':

    options = _import_settings()
    m = pyo.ConcreteModel()

    # Test full model creation with discretization
    print('\nTesting Infinite Horizon Model Setup')
    m = _make_finite_horizon_model(m, options)
    # try:
    # _make_finite_horizon_model(m, options)
    # print('Model Setup and Discretization Successful')
    # except Exception as e:
    #     print(f'Model Setup Failed: {e}')

    if True:  # Toggle to False to disable model display
        with open("model_output.txt", "w") as f:
            m.pprint(ostream=f)

    assert isinstance(m, pyo.ConcreteModel)
