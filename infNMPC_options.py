class Options:
    """
    Simulation settings for ternary distillation two-column system NMPC.

    System: Two distillation columns in series separating A, B, C
    - Column 1: Separates A from B+C
    - Column 2: Separates B from C

    States: 246 total (41 trays × 2 columns × 3 variables per tray)
    MVs: 8 (VB1, LT1, D1, B1, VB2, LT2, D2, B2)
    CVs: 3 (product purities)
    """

    def __init__(self):
        # Simulation control
        self.num_horizons = 25               # Number of MPC steps
        self.nfe_finite = 2                 # Finite elements in finite horizon
        self.ncp_finite = 3                  # Collocation points per FE (finite)
        self.sampling_time = 1               # Time between MPC updates (hours)

        # Infinite horizon settings
        self.infinite_horizon = True         # Use infinite horizon NMPC
        self.nfe_infinite = 1                # Finite elements in infinite horizon
        self.ncp_infinite = 3                # Collocation points per FE (infinite)

        # Solver and model options
        self.tee_flag = True           # Print solver output
        self.endpoint_constraints = False     # Enforce endpoint constraints
        self.custom_objective = True         # Use economic objective
        self.initialize_with_initial_data = True
        self.terminal_cost_riemann = True
        self.remove_collocation = True
        self.initialization_assist = False

        self.input_suppression = False
        self.input_suppression_factor = 1.0e3  # Penalty for input changes

        # Stage cost weights for tracking objective
        # Order: [x1[41,1], x2[41,2], xC, VB1, LT1, D1, B1, VB2, LT2, D2, B2]
        # CVs first, then MVs
        self.stage_cost_weights = [
            1.0e4,  # x1[41,1] - A purity weight
            1.0e4,  # x2[41,2] - B purity weight
            1.0e4,  # xC - C purity weight
            1.0,    # VB1 weight
            1.0,    # LT1 weight
            1.0,    # D1 weight
            1.0,    # B1 weight
            1.0,    # VB2 weight
            1.0,    # LT2 weight
            1.0,    # D2 weight
            1.0,    # B2 weight
        ]

        self.gamma = 0.05                    # Infinite horizon decay parameter
        self.beta = 1.0                      # Terminal cost weight

        # Display/Data Output options
        self.live_plot = True               # Real-time plotting
        self.plot_end = True                 # Plot at end of simulation
        self.save_data = True                # Save simulation data
        self.save_figure = True              # Save plots to file

        # Disturbance options
        self.disturb_flag = False            # Enable process disturbances
        self.disturb_distribution = 'normal'
        self.disturb_seeded = True

        # Slack variable options
        # When includes_slacks=True, the objective penalises slack variables at
        # interior collocation points (those between FE endpoints) in addition
        # to the standard FE-endpoint stage cost. This prevents the degeneracy
        # where slacks at interior points are outside the objective and can grow
        # without bound, corrupting the simulation trajectory.
        self.includes_slacks = True
        self.slack_index = [
            "x1eps", "M1eps", "V1eps", "L1eps", "D1eps", "B1eps", "TC1eps",
            "x2eps", "M2eps", "V2eps", "L2eps", "D2eps", "B2eps", "TC2eps",
            "xCeps",
        ]
        self.slack_rho = 1.0e4    # Penalty weight — keep in sync with rho in custom_objective


def _import_settings():
    """
    Create and return an Options object containing all NMPC settings.

    Returns:
    --------
    Options
        An instance of the Options class with initialized values.
    """
    return Options()
