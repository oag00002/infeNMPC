import pyomo.environ as pyo
from make_model import _make_infinite_horizon_model, _make_finite_horizon_model, _remove_non_collocation_values_infinite, _remove_non_collocation_values_finite
import idaes
from idaes.core.solvers import use_idaes_solver_configuration_defaults
from infNMPC_options import _import_settings
from tqdm import tqdm
from data_save_and_plot import _finalize_live_plot, _setup_live_plot, _update_live_plot, _handle_mpc_results, _save_epsilon, _plot_lyap, _save_lyap_csv
from indexing_tools import _add_time_indexed_expression
import numpy as np
from pyomo.contrib.mpc import ScalarData
from indexing_tools import _get_variable_key_for_data
import time
from math import isclose
from initialization_tools import _assist_initialization_infinite, _assist_initialization_finite
from controller_factory import _make_infinite_horizon_controller, _make_finite_horizon_controller, _make_plant

# Solver Settings
use_idaes_solver_configuration_defaults()
idaes.cfg.ipopt.options.linear_solver = "ma57" # "ma27"
idaes.cfg.ipopt.options.OF_ma57_automatic_scaling = "yes"
idaes.cfg.ipopt.options.max_iter = 6000
idaes.cfg.ipopt.options.halt_on_ampl_error = "yes"
idaes.cfg.ipopt.options.bound_relax_factor = 1e-8
idaes.cfg.ipopt.options.tol = 1e-4
idaes.cfg.ipopt.options.acceptable_tol = 1e-4
# idaes.cfg.ipopt.options.expect_infeasible_problem = "no"
# idaes.cfg.ipopt.options.bound_push = 1e-6
# idaes.cfg.ipopt.options.bound_frac = 1e-6
# idaes.cfg.ipopt.options.mu_strategy = "monotone"
# idaes.cfg.ipopt.options.warm_start_init_point = "yes"
# idaes.cfg.ipopt.options.warm_start_bound_push = 1e-6
# idaes.cfg.ipopt.options.warm_start_mult_bound_push = 1e-6


def _mpc_loop(options):
    """
    Execute the NMPC closed-loop simulation.

    Runs the nonlinear model predictive control loop using either a finite or infinite
    horizon controller. It interacts with the plant model at each sampling interval,
    solves optimization problems, updates plots, and logs performance data.

    Args:
        options (Namespace): Runtime configuration including control horizon type,
            sampling time, plotting options, and solver verbosity.

    Returns:
        None
    """
    if options.infinite_horizon:
        if options.initialization_assist:
            controller = _assist_initialization_infinite(options)
            t0_controller = controller.finite_block.time.first()
            state_vars = controller.finite_block.state_vars
        else:
            controller = _make_infinite_horizon_controller(options)
            t0_controller = controller.finite_block.time.first()
            state_vars = controller.finite_block.state_vars
    else:
        if options.initialization_assist:
            controller = _assist_initialization_finite(options)
            t0_controller = controller.time.first()
            state_vars = controller.state_vars
        else:
            controller = _make_finite_horizon_controller(options)
            t0_controller = controller.time.first()
            state_vars = controller.state_vars

    plant = _make_plant(options)

    idaes.cfg.ipopt.options.warm_start_init_point = "yes"
    idaes.cfg.ipopt.options.warm_start_bound_push = 1e-6
    idaes.cfg.ipopt.options.warm_start_mult_bound_push = 1e-6

    # Get full initial TimeSeriesData
    sim_data = plant.interface.get_data_at_time([plant.time.first()])

    # Zero out non-state values in initial sim_data
    state_var_keys = set()
    for sv in state_vars:
        sv_name = sv.name.split(".")[-1]
        if sv.is_indexed():
            for index in sv.index_set():
                if isinstance(index, tuple) and len(index) > 1:
                    partial_index_str = ",".join(str(i) for i in index[:-1]) + ",*"
                    key = f"{sv_name}[{partial_index_str}]"
                    state_var_keys.add(key)
                else:
                    key = f"{sv_name}[*]"
                    state_var_keys.add(key)
        else:
            key = f"{sv_name}[*]"
            state_var_keys.add(key)

    state_var_keys_str = {str(k) for k in state_var_keys}

    for key in sim_data._data:
        if str(key) not in state_var_keys_str:
            sim_data._data[key][0] = None  # Nullify non-state values

    new_data_time = list(plant.time)[1:]
    solver = pyo.SolverFactory('ipopt')

    # Prepare plotting data arrays
    io_data_array = []
    time_series = []
    solver_time = []
    if options.custom_objective:
        lyap = []

    # Get initial CV values only if they correspond to a differential state
    differential_state_keys = set()
    for sv in state_vars:
        sv_name = sv.name.split(".")[-1]
        if sv.is_indexed():
            for index in sv.index_set():
                if isinstance(index, tuple) and len(index) > 1:
                    partial_index_str = ",".join(str(i) for i in index[:-1]) + ",*"
                    key = f"{sv_name}[{partial_index_str}]"
                    differential_state_keys.add(key)
                else:
                    key = f"{sv_name}[*]"
                    differential_state_keys.add(key)
        else:
            key = f"{sv_name}[*]"
            differential_state_keys.add(key)

    initial_cv_row = []
    for var_name in plant.CV_index:
        key = _get_variable_key_for_data(plant, var_name)
        value = sim_data.get_data_from_key(key)
        if isinstance(value, list):
            initial_cv_row.extend(value)
        else:
            initial_cv_row.append(value)

    initial_mv_row = [None for _ in plant.MV_index]  # No initial MV values
    io_data_array.append(initial_cv_row + initial_mv_row)
    time_series.append(0.0)

    fig, axes = None, None
    if options.live_plot:
        fig, axes = _setup_live_plot(plant)

    loop_iter = tqdm(range(options.num_horizons), desc="Running MPC")

    if options.infinite_horizon:
        terminal_cost_prev = 1
        first_stage_cost_prev = 1

    for i in loop_iter:

        simulation_time = (i + 1) * options.sampling_time

        start_time = time.perf_counter()
        solver.solve(controller, tee=options.tee_flag)
        end_time = time.perf_counter()

        # if i == 1:  # Toggle to False to disable model display
        #     with open("model_output.txt", "w") as f:
        #         controller.pprint(ostream=f)

        solver_time.append(end_time - start_time)

        if options.infinite_horizon:

            # Get all time points and finite element boundaries
            t_points = controller.finite_block.time
            fe_points = t_points.get_finite_elements()

            if options.endpoint_constraints:
                terminal_cost_now = pyo.value(
                    controller.infinite_block.phi[controller.infinite_block.time.last()]
                    * options.beta / options.sampling_time
                )
            else:
                terminal_cost_now = pyo.value(
                    controller.infinite_block.phi[controller.infinite_block.time.prev(controller.infinite_block.time.last())]
                )

            c = options.stage_cost_weights

            penultimate_stage_cost_now = sum(
                pyo.value(
                    c[i]
                    * (
                        _add_time_indexed_expression(
                            controller.finite_block,
                            var_name,
                            fe_points[-1]
                        ) - controller.finite_block.steady_state_values[var_name]
                    )**2
                )
                for i, var_name in enumerate(controller.finite_block.stage_cost_index)
            )

            # Optional debug printing
            print("")
            print("Terminal cost (prev):", terminal_cost_prev)
            print("Terminal cost (now):", terminal_cost_now)
            print("First stage cost (prev):", first_stage_cost_prev)
            print("Penultimate stage cost (now):", penultimate_stage_cost_now)
            print("Lyapunov Function Value:", pyo.value(controller.lyapunov))
            lyap.append(pyo.value(controller.lyapunov))
            if i == 0:
                lyap.append(pyo.value(controller.lyapunov))

            dlyap = lyap[i+1] - lyap[i]

            LHS = -(
                terminal_cost_prev - terminal_cost_now - options.beta * penultimate_stage_cost_now
            ) / first_stage_cost_prev

            # if simulation_time != options.sampling_time:
            #     assert LHS < 1, f"No ε ∈ [0,1) satisfies LHS ≤ ε; got LHS = {LHS}"
            print("")
            print(f"Min value of epsilon: {LHS}")
            print(f"Change in Lyapunov Function Value: {dlyap}")
            print("")

            custom_obj = terminal_cost_prev + first_stage_cost_prev - terminal_cost_now - penultimate_stage_cost_now

            custom_obj *= -1

            if options.save_data:
                _save_epsilon(i, custom_obj, options)
        

        if options.infinite_horizon:
            if options.endpoint_constraints:
                terminal_cost_prev = pyo.value(
                    controller.infinite_block.phi[controller.infinite_block.time.last()]
                    * options.beta / options.sampling_time
                )
            else:
                terminal_cost_prev = pyo.value(
                    controller.infinite_block.phi[controller.infinite_block.time.prev(controller.infinite_block.time.last())]
                    * options.beta / options.sampling_time
                )

            first_stage_cost_prev = sum(
                pyo.value(
                    c[i]
                    * (
                        _add_time_indexed_expression(
                            controller.finite_block,
                            var_name,
                            fe_points[1]
                        ) - controller.finite_block.steady_state_values[var_name]
                    )**2
                )
                for i, var_name in enumerate(controller.finite_block.stage_cost_index)
            )

        ts_data = controller.interface.get_data_at_time(options.sampling_time)

        if options.infinite_horizon:
            input_data = ts_data.extract_variables(
                [getattr(controller.finite_block, var_name) for var_name in controller.finite_block.MV_index],
                context=controller.finite_block
            )
        else:
            input_data = ts_data.extract_variables(
                [getattr(controller, var_name) for var_name in controller.MV_index]
            )

        plant.interface.load_data(input_data, time_points=new_data_time)
        solver.solve(plant, tee=options.tee_flag)

        # Get CV and MV values at the sampling time
        full_data = plant.interface.get_data_at_time(options.sampling_time)
        cv_row = []
        for var_name in plant.CV_index:
            val = full_data.get_data_from_key(_get_variable_key_for_data(plant, var_name))
            if isinstance(val, list):
                cv_row.extend(val)
            else:
                cv_row.append(val)

        mv_row = []
        for var_name in plant.MV_index:
            val = full_data.get_data_from_key(_get_variable_key_for_data(plant, var_name))
            if isinstance(val, list):
                mv_row.extend(val)
            else:
                mv_row.append(val)

        io_row = cv_row + mv_row
        io_data_array.append(io_row)
        time_series.append(simulation_time)

        if options.live_plot:
            _update_live_plot(fig, axes, time_series, io_data_array, plant)

        model_data = plant.interface.get_data_at_time(new_data_time)
        model_data.shift_time_points(simulation_time - plant.time.first() - options.sampling_time)
        sim_data.concatenate(model_data)

        # Only load state variables into tf_data
        full_tf_data = plant.interface.get_data_at_time(options.sampling_time)
        tf_data = ScalarData(data={})
        for state_var in state_vars:
            sv_name = state_var.name.split(".")[-1]
            if state_var.is_indexed():
                for index in state_var.index_set():
                    if isinstance(index, tuple) and len(index) > 1:
                        partial_index_str = ",".join(str(i) for i in index[:-1]) + ",*"
                        key = f"{sv_name}[{partial_index_str}]"
                        tf_data._data[key] = full_tf_data.get_data_from_key(key)
                    else:
                        key = f"{sv_name}[*]"
                        tf_data._data[key] = full_tf_data.get_data_from_key(key)
            else:
                key = f"{sv_name}[*]"
                tf_data.data[key] = full_tf_data.get_data_from_key(key)

        plant.interface.load_data(tf_data)

        controller.interface.shift_values_by_time(options.sampling_time)
        controller.interface.load_data(tf_data, time_points=t0_controller)
        
        if options.remove_collocation:
            if options.infinite_horizon:
                controller = _remove_non_collocation_values_infinite(controller)
            else:
                controller = _remove_non_collocation_values_finite(controller)

    if options.custom_objective:
        if options.plot_end:
            _plot_lyap(lyap, options)
        if options.save_data:
            _save_lyap_csv(lyap, options)

    if options.live_plot and fig is not None:
        _finalize_live_plot(fig)

    _handle_mpc_results(sim_data, time_series, io_data_array, plant, solver_time, options)


if __name__ == "__main__":
    options = _import_settings()
    _mpc_loop(options)