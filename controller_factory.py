import pyomo.environ as pyo
from make_model import _make_infinite_horizon_model, _make_finite_horizon_model
from indexing_tools import _add_time_indexed_expression

def _make_infinite_horizon_controller(options, data=None):
    """
    Create and solve a Pyomo model for infinite-horizon NMPC.

    This function builds an infinite-horizon nonlinear model predictive controller
    using a model defined in `make_model.py`, constructs a stage and terminal cost
    objective, and solves the model using IPOPT.

    Args:
        options (Namespace): Configuration options for the controller, including
            weights, solver settings, and model parameters.

    Returns:
        pyo.ConcreteModel: The solved infinite-horizon Pyomo model instance.
    """
    m = pyo.ConcreteModel()
    m = _make_infinite_horizon_model(m, options)

    if data is not None:
        # Load initial data into the model
        m.interface.load_data(data, time_points=list(m.time))

    print('Generating Infinite Horizon Controller')

    m.finite_block.stage_cost_index = pyo.Set(initialize=lambda m: list(m.CV_index) + list(m.MV_index))

    # Define the objective function
    if options.custom_objective:
        from model_equations import custom_objective

        # get the time-dependent cost function
        cost_fn = custom_objective(m.finite_block, options)

        def objective_rule(m):
            finite_elements = m.finite_block.time.get_finite_elements()
            # Sum the custom stage cost across the finite time horizon
            stage_cost = sum(cost_fn(m.finite_block, t) - m.finite_block.ss_obj_value for t in finite_elements if t != m.finite_block.time.first())

            if options.endpoint_constraints:
                terminal_cost = m.infinite_block.phi[m.infinite_block.time.last()] * options.beta / options.sampling_time
            else:
                # Add terminal cost (same as non-custom case)
                terminal_cost = m.infinite_block.phi[m.infinite_block.time.prev(m.infinite_block.time.last())] * options.beta / options.sampling_time

            return stage_cost + terminal_cost

    elif options.terminal_cost_riemann:
        def objective_rule(m):
            c = options.stage_cost_weights
            finite_elements = m.finite_block.time.get_finite_elements()
            stage_cost = sum(sum(
                c[i] * (_add_time_indexed_expression(m.finite_block, var_name, t) - m.finite_block.steady_state_values[var_name])**2
                for i, var_name in enumerate(m.finite_block.stage_cost_index)
            ) for t in finite_elements if t != m.finite_block.time.first())

            tau_list = list(m.infinite_block.time.data())
            tau_prev = 0
            obj_expression_infinite = 0

            for tau in tau_list:
                if m.infinite_block.time.first() < tau < m.infinite_block.time.last():
                    tracking_cost = sum(
                        c[i] * (_add_time_indexed_expression(m.infinite_block, var_name, tau) - m.infinite_block.steady_state_values[var_name])**2
                        for i, var_name in enumerate(m.infinite_block.stage_cost_index)
                    )
                    obj_expression_infinite += tracking_cost * (tau - tau_prev) / (options.gamma / options.sampling_time * (1 - tau**2))
                    tau_prev = tau

            terminal_cost = obj_expression_infinite * options.beta / options.sampling_time

            return stage_cost + terminal_cost
    else:

        def objective_rule(m):
            c = options.stage_cost_weights
            finite_elements = m.finite_block.time.get_finite_elements()
            stage_cost = sum(sum(
                c[i] * (_add_time_indexed_expression(m.finite_block, var_name, t) - m.finite_block.steady_state_values[var_name])**2
                for i, var_name in enumerate(m.finite_block.stage_cost_index)
            ) for t in finite_elements if t != m.finite_block.time.first())

            if options.input_suppression:
                for t in m.finite_block.time:
                    if t > m.finite_block.time.first():
                        t_prev = m.finite_block.time.prev(t)
                        for mv in m.finite_block.MV_index:
                            mv_expr = _add_time_indexed_expression(m.finite_block, mv, t) - _add_time_indexed_expression(m.finite_block, mv, t_prev)
                            stage_cost += options.input_suppression_factor * mv_expr**2
                
                for t in m.infinite_block.time:
                    if t > m.infinite_block.time.first():
                        t_prev = m.infinite_block.time.prev(t)
                        for mv in m.infinite_block.MV_index:
                            mv_expr = _add_time_indexed_expression(m.infinite_block, mv, t) - _add_time_indexed_expression(m.infinite_block, mv, t_prev)
                            stage_cost += options.input_suppression_factor * mv_expr**2

            if options.endpoint_constraints:
                terminal_cost = m.infinite_block.phi[m.infinite_block.time.last()] * options.beta / options.sampling_time
            else:
                # Add terminal cost (same as non-custom case)
                terminal_cost = m.infinite_block.phi[m.infinite_block.time.prev(m.infinite_block.time.last())] * options.beta / options.sampling_time

            return stage_cost + terminal_cost

    m.objective = pyo.Objective(rule=objective_rule, sense=pyo.minimize)

    print('Infinite Horizon Controller Initial Solve')

    solver = pyo.SolverFactory('ipopt')
    solver.solve(m, tee=options.tee_flag)

    if True:  # Toggle to False to disable model display
        with open("model_output.txt", "w") as f:
            m.pprint(ostream=f)

    return m


def _make_finite_horizon_controller(options, data=None):
    """
    Create and solve a Pyomo model for finite-horizon NMPC.

    Builds a finite-horizon nonlinear MPC model, optionally using a custom objective 
    or a standard quadratic tracking objective. Solves the resulting optimization 
    problem with IPOPT.

    Args:
        options (Namespace): Controller configuration including weights, solver options,
            and custom objective flag.

    Returns:
        pyo.ConcreteModel: The solved finite-horizon Pyomo model instance.
    """
    m = pyo.ConcreteModel()
    m = _make_finite_horizon_model(m, options)

    if data is not None:
        print("Loading initial data into finite horizon controller")
        for i in range(len(data)):
            m.interface.load_data(data[i], time_points=list(m.time)[i])

    if True:  # Toggle to False to disable model display
        with open("model_output.txt", "w") as f:
            m.pprint(ostream=f)

    print('Generating Finite Horizon Controller')

    m.stage_cost_index = pyo.Set(initialize=lambda m: list(m.CV_index) + list(m.MV_index))

    # Define the objective function
    if options.custom_objective:
        from model_equations import custom_objective

        # get the time-dependent cost function
        cost_fn = custom_objective(m, options)

        def objective_rule(m):
            # Sum the custom stage cost across the finite time horizon
            finite_elements = m.time.get_finite_elements()
            stage_cost = sum(cost_fn(m, t) - m.ss_obj_value for t in finite_elements if t != m.time.first())

            return stage_cost

    else:
        def objective_rule(m):
            c = options.stage_cost_weights
            finite_elements = m.time.get_finite_elements()
            
            stage_cost = sum(sum(
                c[i] * (_add_time_indexed_expression(m, var_name, t) - m.steady_state_values[var_name])**2
                for i, var_name in enumerate(m.stage_cost_index)
            ) for t in finite_elements if t != m.time.first()) #  for t in m.time if t != m.time.first()) #

            if options.input_suppression:
                for t in m.time:
                    if t > m.time.first():
                        t_prev = m.time.prev(t)
                        for mv in m.MV_index:
                            mv_expr = _add_time_indexed_expression(m, mv, t) - _add_time_indexed_expression(m, mv, t_prev)
                            stage_cost += options.input_suppression_factor * mv_expr**2
            
            return stage_cost

    m.objective = pyo.Objective(rule=objective_rule, sense=pyo.minimize)

    print('Finite Horizon Controller Initial Solve')

    solver = pyo.SolverFactory('ipopt')
    solver.solve(m, tee=options.tee_flag)

    if True:  # Toggle to False to disable model display
        with open("model_output.txt", "w") as f:
            m.pprint(ostream=f)

    return m
    ######################################################
    # m = pyo.ConcreteModel()
    # m = _make_finite_horizon_model(m, options)

    # m.stage_cost_index = pyo.Set(initialize=lambda m: list(m.CV_index) + list(m.MV_index))

    # if False:
    #     def objective_rule(m_controller):
    #         # -------- Finite Horizon Cost (negative as in GAMS) --------
    #         obj_expression_finite = sum(
    #             -(
    #                 m_controller.Fa0[t] * (2 * m_controller.Cb[t] - 0.5) - 6
    #             ) * (t - m_controller.time.prev(t))
    #             for t in m_controller.time
    #             if t != m_controller.time.first()
    #         )

    #         return obj_expression_finite
    # else:
    #     # Define the objective function
    #     def objective_rule(m):
    #         c = options.stage_cost_weights
    #         stage_cost = sum(sum(
    #             c[i] * (_add_time_indexed_expression(m, var_name, t) - m.steady_state_values[var_name])**2
    #             for i, var_name in enumerate(m.stage_cost_index)
    #         ) for t in m.time)

    #         return stage_cost

    # m.objective = pyo.Objective(rule=objective_rule, sense=pyo.minimize)

    # solver = pyo.SolverFactory('ipopt')
    # solver.solve(m, tee=options.tee_flag)

    # return m