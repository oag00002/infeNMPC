from infNMPC_options import _import_settings
from controller_factory import _make_infinite_horizon_controller, _make_finite_horizon_controller
from pyomo.contrib.mpc.data.scalar_data import ScalarData
from pyomo.contrib.mpc.data.series_data import TimeSeriesData

def _assist_initialization_infinite(options):
    print("Using initialization assist for infinite horizon controller.")

    # Clone options and start with a larger sampling time
    modified_options = _import_settings()
    modified_options.__dict__.update(vars(options))
    modified_options.sampling_time = options.initialization_assist_sampling_time_start

    # Create initial controller and extract data
    controller = _make_infinite_horizon_controller(modified_options)
    interface = controller.interface
    time_points = list(interface.model.finite_block.time) + list(interface.model.infinite_block.time)
    print("Initial time points:", time_points)
    previous_time = list(interface.model.time)
    solution_data_list = [interface.get_data_at_time(t) for t in previous_time]

    # Iteratively reduce sampling time
    while modified_options.sampling_time * 0.9 > options.sampling_time:
        modified_options.sampling_time *= 0.9
        modified_options.sampling_time = round(modified_options.sampling_time, 8)
        print(f"  Trying sampling time: {modified_options.sampling_time}")

        controller = _make_infinite_horizon_controller(modified_options)
        interface = controller.interface
        current_time = list(interface.model.finite_block.time) + list(interface.model.infinite_block.time)
        print("Updated time points:", current_time)
        # Repack solution data as TimeSeriesData aligned by index
        new_data = TimeSeriesData(
            {cuid: [d[cuid] for d in solution_data_list] for cuid in solution_data_list[0]},
            time=current_time,
            time_set=interface.time,
        )
        interface.load_data(new_data, time_points=current_time)
        solution_data_list = [interface.get_data_at_time(t) for t in interface.model.time]

    # Final solve at target sampling time
    print(f"Final solve at target sampling time: {options.sampling_time}")
    controller = _make_infinite_horizon_controller(options)
    interface = controller.interface
    final_time = list(interface.model.finite_block.time) + list(interface.model.infinite_block.time)
    final_data = TimeSeriesData(
        {cuid: [d[cuid] for d in solution_data_list] for cuid in solution_data_list[0]},
        time=final_time,
        time_set=interface.time,
    )
    interface.load_data(final_data, time_points=final_time)
    options.final_controller = controller

def _assist_initialization_finite(options):
    print("Using initialization assist for finite horizon controller.")

    # Clone options and start with a larger sampling time
    modified_options = _import_settings()
    modified_options.__dict__.update(vars(options))
    modified_options.sampling_time = options.initialization_assist_sampling_time_start

    # Create initial controller and extract data
    m = _make_finite_horizon_controller(modified_options)
    previous_time = list(m.time)
    print("Initial time points:", previous_time)
    solution_data_list = [m.interface.get_data_at_time(t) for t in previous_time]

    # Iteratively reduce sampling time
    while modified_options.sampling_time * 0.9 > options.sampling_time:
        modified_options.sampling_time *= 0.9
        modified_options.sampling_time = round(modified_options.sampling_time, 8)
        print(f"  Trying sampling time: {modified_options.sampling_time}")

        m = _make_finite_horizon_controller(modified_options, solution_data_list)
        current_time = list(m.time)
        print("Updated time points:", current_time)
        solution_data_list = [m.interface.get_data_at_time(t) for t in current_time]

    # Final solve at target sampling time
    print(f"Final solve at target sampling time: {options.sampling_time}")
    m = _make_finite_horizon_controller(options, solution_data_list)
    return m