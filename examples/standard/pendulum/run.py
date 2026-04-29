"""
Run script for the cart-pole (inverted pendulum) example.

System:  Cart-pole
States:  x (m), x_dot (m/s), theta (rad), theta_dot (rad/s)
MV:      F — horizontal force on cart (N)
Objective: quadratic regulation of theta and x to zero

Usage
-----
From the repository root:

    python examples/pendulum/run.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

import model  # noqa: E402

from infeNMPC import Options, mpc_loop  # noqa: E402

options = Options.for_model_module(
    model,
    # ---- Simulation control ----
    num_horizons=150,
    sampling_time=0.2,
    # ---- Finite-horizon discretization ----
    nfe_finite=5,
    ncp_finite=5,
    # ---- Infinite-horizon discretization ----
    nfe_infinite=6,
    ncp_infinite=5,
    # ---- Controller settings ----
    infinite_horizon=True,
    terminal_constraint_type='hard',
    objective='tracking',
    initialize_with_initial_data=True,
    # ---- Cost function ----
    stage_cost_weights=[1, 1, 1],
    gamma=0.0075,
    beta=1.2,
    # ---- Output ----
    tee_flag=False,
    save_data=True,
    save_figure=True,
)

if __name__ == '__main__':
    mpc_loop(options)
