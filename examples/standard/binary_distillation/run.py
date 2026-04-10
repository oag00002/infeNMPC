"""
Run script for the binary distillation column example (setpoint tracking).

System:  42-tray Methanol/n-Propanol distillation column
States:  M[tray] (holdup), x[tray] (mole fraction)
MVs:     Qr — reboiler heat duty (MJ/h), Rec — reflux ratio
CVs:     x[1] — bottoms composition, T[29] — mid-column temperature
Objective: quadratic tracking to steady-state setpoints

Usage
-----
From the repository root:

    python examples/binary_distillation/run.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

import model  # noqa: E402

from infeNMPC import Options, mpc_loop  # noqa: E402

options = Options.for_model_module(
    model,
    # ---- Simulation control ----
    num_horizons=25,
    sampling_time=10,
    # ---- Finite-horizon discretization ----
    nfe_finite=25,
    ncp_finite=3,
    # ---- Infinite-horizon discretization ----
    nfe_infinite=3,
    ncp_infinite=3,
    # ---- Controller settings ----
    infinite_horizon=False,
    terminal_constraint_type='hard',
    custom_objective=False,
    # ---- Cost function ----
    stage_cost_weights=[1e5, 1e3, 0, 0],
    gamma=0.05,
    beta=1.0,
    # ---- Output ----
    tee_flag=True,
    save_data=True,
    save_figure=True,
)

if __name__ == '__main__':
    mpc_loop(options)
