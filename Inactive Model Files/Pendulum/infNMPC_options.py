class Options:
    """
    A class to store and manage simulation settings for infinite and finite horizon NMPC.

    Attributes:
    -----------
    num_horizons : int
        Number of simulation steps to run.
    finite_horizon : float
        Time length of the finite horizon (in time units, e.g., hours or seconds).
    ncp_finite : int
        Number of collocation points per finite-horizon finite element.
    sampling_time : float
        Time interval between successive measurements or control updates.
    infinite_horizon : bool
        If True, infinite horizon MPC will be used. If False, finite horizon MPC is used.
    nfe_infinite : int
        Number of finite elements for the infinite horizon approximation.
    ncp_infinite : int
        Number of collocation points per infinite-horizon finite element.
    tee_flag : bool
        If True, solver output will be printed to the console.
    finite_terminal_con : bool
        If True, a terminal state constraint is applied in the finite horizon case.
    terminal_cost_riemann : bool
        If True, a Riemann sum is used to approximate the terminal cost.
    """

    def __init__(self):
        # Simulation control
        self.num_horizons = 150
        self.nfe_finite = 5
        self.ncp_finite = 5
        self.sampling_time = 0.2

        # Infinite horizon settings
        self.infinite_horizon = True
        self.nfe_infinite = 6
        self.ncp_infinite = 5

        # Solver and model options
        self.tee_flag = True
        self.endpoint_constraints = True
        self.custom_objective = False
        self.initialize_with_initial_data = True
        self.terminal_cost_riemann = False

        self.stage_cost_weights = [1, 1, 1]
        self.gamma = 0.0075
        self.beta = 1.2

        # Display/Data Output options
        self.live_plot = True
        self.plot_end = False
        self.save_data = True
        self.save_figure = True

        # Disturbance options
        self.disturb_flag = False
        self.disturb_distribution = 'normal'
        self.disturb_seeded = True


def _import_settings():
    """
    Create and return an Options object containing all NMPC settings.

    Returns:
    --------
    Options
        An instance of the Options class with initialized values.
    """
    return Options()