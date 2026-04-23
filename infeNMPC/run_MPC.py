import os
import pyomo.environ as pyo
from tqdm import tqdm
from pyomo.contrib.mpc import ScalarData

from .make_model import _ipopt_solver
from .infNMPC_options import Options, _import_settings
from .plant import Plant
from .controllers import InfiniteHorizonController, FiniteHorizonController
from .tools.initialization_tools import _assist_initialization_infinite, _assist_initialization_finite
from .data_save_and_plot import (
    _handle_mpc_results,
    _create_run_folder,
    _save_io_csv,
)
from .tools.indexing_tools import _get_variable_key_for_data, _add_time_indexed_expression
from .tools.debug_tools import _report_constraint_violations


def mpc_loop(options: Options):
    """
    Execute the NMPC closed-loop simulation.

    Runs the nonlinear model predictive control loop using either a finite or
    infinite horizon controller.  Interacts with the plant model at each
    sampling interval, solves optimization problems, updates plots, and logs
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
        t_last_controller = controller.finite_block.time.last()
        state_vars = controller.finite_block.state_vars
        resolved_gamma = controller.infinite_block.gamma
        _ss_ref_block = controller.finite_block
    else:
        if options.initialization_assist:
            controller = _assist_initialization_finite(options)
        else:
            controller = FiniteHorizonController(options)
        t0_controller = controller.time.first()
        t_last_controller = controller.time.last()
        state_vars = controller.state_vars
        resolved_gamma = None
        _ss_ref_block = controller._model

    # Steady-state warm data for the horizon endpoint after each shift.
    # After shift_values_by_time, the last time point gets the old endpoint
    # value repeated (there is no t+h data beyond the horizon).  For a
    # Lyapunov-constrained controller this doubled endpoint inflates phi_track
    # above V_prev, giving IPOPT a large initial infeasibility.  Loading the
    # steady-state values at the last time point resets the guess for the
    # "extrapolated" tail to the economic optimum, keeping phi_track at the
    # warm-started point below V_prev.
    _ss_warm_data = ScalarData({
        _get_variable_key_for_data(_ss_ref_block, var): _ss_ref_block.steady_state_values[var]
        for var in list(_ss_ref_block.CV_index) + list(_ss_ref_block.MV_index)
    })

    # Create the timestamped results folder now so safe_run writes go there too.
    run_folder = _create_run_folder(options, resolved_gamma=resolved_gamma)

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

    safe_run_folder = run_folder if options.safe_run else None

    loop_iter = tqdm(range(options.num_horizons), desc="Running MPC")

    for i in loop_iter:

        if options.debug_flag and i == 0:
            _dbg_dir = os.path.join("research_files", "shifting_behavior")
            os.makedirs(_dbg_dir, exist_ok=True)
            with open(os.path.join(_dbg_dir, "plant_model_init.txt"), "w") as f:
                plant._model.pprint(ostream=f)
            with open(os.path.join(_dbg_dir, "controller_model_pre_solve.txt"), "w") as f:
                controller.pprint(ostream=f)

        simulation_time = (i + 1) * options.sampling_time

        try:
            controller.solve()
        except RuntimeError:
            if options.debug_flag:
                _report_constraint_violations(
                    controller._model, label=f"controller failure at iteration {i}"
                )
            raise
        cpu_time.append(controller.last_solve_time)

        if options.debug_flag and i == 0:
            _shift_obs_dir = os.path.join("research_files", "shifting_behavior")
            os.makedirs(_shift_obs_dir, exist_ok=True)
            with open(os.path.join(_shift_obs_dir, f"iter{i:03d}_post_solve.txt"), "w") as f:
                controller.pprint(ostream=f)

        # Update Lyapunov stability constraint parameters for next iteration
        if options.lyap_flag:
            # V_prev / first_stage_cost_prev live on controller._model in both cases.
            # For infinite horizon: V = finite Riemann sum + infinite ODE integral.
            # For finite horizon:   V = finite Riemann sum only.
            lyap_block = controller
            if options.infinite_horizon:
                stage_block = controller.finite_block
                lyap_beta = getattr(options, 'lyap_beta', 1.0)
                V_current = (
                    pyo.value(controller.finite_block.phi_track)
                    + lyap_beta * pyo.value(
                        controller.infinite_block.phi_track[
                            controller.infinite_block.time.last()
                        ]
                    )
                )
            else:
                stage_block = controller
                V_current = pyo.value(controller.phi_track)
            track_list = list(stage_block.CV_index) + list(stage_block.MV_index)
            c_raw = options.stage_cost_weights or []
            c_track = c_raw[:len(track_list)] if len(c_raw) >= len(track_list) else [1.0] * len(track_list)
            t_first_fe = next(
                t for t in stage_block.time.get_finite_elements()
                if t > stage_block.time.first()
            )
            first_stage_cost = sum(
                c_track[j] * (
                    pyo.value(_add_time_indexed_expression(stage_block, var, t_first_fe))
                    - stage_block.steady_state_values[var]
                ) ** 2
                for j, var in enumerate(track_list)
            )
            lyap_block.V_prev.set_value(V_current)
            lyap_block.first_stage_cost_prev.set_value(first_stage_cost)

        # Warm-start the plant with the controller's full first-FE solution at every
        # collocation point.  Since the plant time domain is identical to the
        # controller's first FE (same ncp, same sampling_time), controller.interface
        # has values for exactly the time points in new_data_time.  Loading all
        # variables (state, CV, MV, slack) means:
        #   (a) the plant starts at the controller's solution → 1-2 IPOPT iterations
        #       in the no-disturbance case.
        #   (b) MV and slack variables receive the controller-optimised values via
        #       load_data.  Both are fixed in the plant (MVs from Plant.__init__,
        #       slacks fixed below), so they must be temporarily unfixed to allow
        #       load_data to update them, then re-fixed after loading.
        plant_warmstart = controller.interface.get_data_at_time(new_data_time)

        # Temporarily unfix MVs and slacks so load_data can propagate
        # controller values.  load_data skips fixed variables, so both must
        # be free before the call.  MVs are fixed from Plant.__init__; slacks
        # are fixed at the end of each iteration (including the first, where
        # Plant.__init__ leaves them free).  Re-fix both after loading.
        for _mv_name in plant._model.MV_index:
            _mv = getattr(plant._model, _mv_name, None)
            if _mv is not None:
                _mv.unfix()
        if hasattr(plant._model, 'slack_index'):
            for _sv_name in plant._model.slack_index:
                _sv = getattr(plant._model, _sv_name, None)
                if _sv is not None and isinstance(_sv, pyo.Var):
                    _sv.unfix()

        plant.interface.load_data(plant_warmstart, time_points=new_data_time)

        # Re-fix MVs and slacks at their newly loaded values.
        for _mv_name in plant._model.MV_index:
            _mv = getattr(plant._model, _mv_name, None)
            if _mv is not None:
                _mv.fix()
        if hasattr(plant._model, 'slack_index'):
            for _sv_name in plant._model.slack_index:
                _sv = getattr(plant._model, _sv_name, None)
                if _sv is not None and isinstance(_sv, pyo.Var):
                    _sv.fix()

        if options.debug_flag and i <= 1:
            _dbg_dir = os.path.join("research_files", "shifting_behavior")
            os.makedirs(_dbg_dir, exist_ok=True)
            with open(os.path.join(_dbg_dir, f"iter{i:03d}_controller_post_solve.txt"), "w") as f:
                controller.pprint(ostream=f)
            with open(os.path.join(_dbg_dir, f"iter{i:03d}_plant_pre_solve.txt"), "w") as f:
                plant._model.pprint(ostream=f)

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

        if options.safe_run:
            _save_io_csv(time_series, io_data_array, plant, safe_run_folder)

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

        # Load state-var IC into plant at t=0.
        # All state vars at t=0 are fixed (no algebraic CV coupling constraints
        # exist there after _finite_block_gen deletes equations_write t=0 entries),
        # so load_data skips them.  Temporarily unfix, load, then refix.
        plant.interface.load_data(tf_data)  # warm-start at t>0 (t=0 skipped, fixed)
        _t0_plant = plant._model.time.first()
        for _sv in plant._model.state_vars:
            for _idx in _sv.index_set():
                if (_idx[-1] if isinstance(_idx, tuple) else _idx) == _t0_plant:
                    _sv[_idx].unfix()
        plant.interface.load_data(tf_data, time_points=_t0_plant)
        for _sv in plant._model.state_vars:
            for _idx in _sv.index_set():
                if (_idx[-1] if isinstance(_idx, tuple) else _idx) == _t0_plant:
                    _sv[_idx].fix()

        # Shift controller warm-start and load exact plant endpoint as new IC.
        # shift_values_by_time updates all vars (including fixed) via set_value,
        # but load_data then skips fixed state vars.  Unfix → load → refix ensures
        # the controller IC exactly matches the plant endpoint rather than the
        # controller's own shifted endpoint.
        controller.interface.shift_values_by_time(options.sampling_time)
        for _sv in state_vars:
            for _idx in _sv.index_set():
                if (_idx[-1] if isinstance(_idx, tuple) else _idx) == t0_controller:
                    _sv[_idx].unfix()
        controller.interface.load_data(tf_data, time_points=t0_controller)
        for _sv in state_vars:
            for _idx in _sv.index_set():
                if (_idx[-1] if isinstance(_idx, tuple) else _idx) == t0_controller:
                    _sv[_idx].fix()
        controller.interface.load_data(_ss_warm_data, time_points=t_last_controller)

        if options.debug_flag and i == 0:
            with open(os.path.join(_shift_obs_dir, f"iter{i:03d}_post_shift.txt"), "w") as f:
                controller.pprint(ostream=f)

    _handle_mpc_results(sim_data, time_series, io_data_array, plant, cpu_time, options,
                        folder_path=run_folder)


if __name__ == "__main__":
    options = _import_settings()
    mpc_loop(options)
