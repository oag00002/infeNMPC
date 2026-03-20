import pyomo.environ as pyo
from tqdm import tqdm
from pyomo.contrib.mpc import ScalarData

from .make_model import _ipopt_solver
from .infNMPC_options import Options, _import_settings
from .plant import Plant
from .controllers import InfiniteHorizonController, FiniteHorizonController
from .initialization_tools import _assist_initialization_infinite, _assist_initialization_finite
from .data_save_and_plot import (
    _finalize_live_plot, _setup_live_plot, _update_live_plot,
    _handle_mpc_results,
)
from .indexing_tools import _get_variable_key_for_data


def mpc_loop(options: Options):
    """
    Execute the NMPC closed-loop simulation.

    Runs the nonlinear model predictive control loop using either a finite or
    infinite horizon controller.  Interacts with the plant model at each
    sampling interval, solves optimisation problems, updates plots, and logs
    performance data.

    Args:
        options (Options): Runtime configuration including control horizon type,
            sampling time, plotting options, and solver verbosity.
    """
    if options.infinite_horizon:
        if options.initialization_assist:
            controller = _assist_initialization_infinite(options)
        else:
            controller = InfiniteHorizonController(options)
        t0_controller = controller.finite_block.time.first()
        state_vars = controller.finite_block.state_vars
    else:
        if options.initialization_assist:
            controller = _assist_initialization_finite(options)
        else:
            controller = FiniteHorizonController(options)
        t0_controller = controller.time.first()
        state_vars = controller.state_vars

    plant = Plant(options)

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
            sim_data._data[key][0] = None

    new_data_time = list(plant.time)[1:]

    # Prepare plotting data arrays
    io_data_array = []
    time_series = []
    cpu_time = []

    # Collect initial CV values
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

    initial_mv_row = [None for _ in plant.MV_index]
    io_data_array.append(initial_cv_row + initial_mv_row)
    time_series.append(0.0)

    fig, axes = None, None
    if options.live_plot:
        fig, axes = _setup_live_plot(plant)

    loop_iter = tqdm(range(options.num_horizons), desc="Running MPC")

    for i in loop_iter:

        simulation_time = (i + 1) * options.sampling_time

        controller.solve()
        cpu_time.append(controller.last_solve_time)

        ts_data = controller.interface.get_data_at_time(options.sampling_time)

        if options.infinite_horizon:
            input_data = ts_data.extract_variables(
                [getattr(controller.finite_block, var_name)
                 for var_name in controller.finite_block.MV_index],
                context=controller.finite_block,
            )
        else:
            input_data = ts_data.extract_variables(
                [getattr(controller, var_name) for var_name in controller.MV_index]
            )

        plant.interface.load_data(input_data, time_points=new_data_time)
        plant.solve()

        # Collect CV and MV values at the sampling time
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

        io_data_array.append(cv_row + mv_row)
        time_series.append(simulation_time)

        if options.live_plot:
            _update_live_plot(fig, axes, time_series, io_data_array, plant)

        model_data = plant.interface.get_data_at_time(new_data_time)
        model_data.shift_time_points(
            simulation_time - plant.time.first() - options.sampling_time
        )
        sim_data.concatenate(model_data)

        # Load only state variables for the next initial condition
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

    if options.live_plot and fig is not None:
        _finalize_live_plot(fig)

    _handle_mpc_results(sim_data, time_series, io_data_array, plant, cpu_time, options)


if __name__ == "__main__":
    options = _import_settings()
    mpc_loop(options)
