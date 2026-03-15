"""
Progressive sampling-time warm-start utilities for robust controller initialisation.
"""
from .controllers import InfiniteHorizonController, FiniteHorizonController
from .infNMPC_options import Options
from pyomo.contrib.mpc.data.series_data import TimeSeriesData


def _assist_initialization_infinite(options: Options) -> InfiniteHorizonController:
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


def _assist_initialization_finite(options: Options) -> FiniteHorizonController:
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
