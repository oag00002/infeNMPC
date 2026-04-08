"""
Run script for the eNMPC binary distillation column example (economic objective).

System:  42-tray Methanol/n-Propanol distillation column
States:  M[tray] (holdup), x[tray] (mole fraction)
MVs:     Qr — reboiler heat duty (MJ/h), Rec — reflux ratio
CVs:     x[1] — bottoms composition, T[29] — mid-column temperature
Objective: maximise bottoms purity while penalising reboiler and condenser duty
    stage_cost = -(100*(1 - x[1,t]) - 10*Qr[t] + Qc[t])

Usage
-----
From the repository root:

    python examples/enmpc_binary_distillation/run.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

import binary_distillation_model  # noqa: E402

from infeNMPC import Options, mpc_loop  # noqa: E402

options = Options.for_model_module(
    binary_distillation_model,
    # ---- Simulation control ----
    num_horizons=25,
    sampling_time=10,
    # ---- Finite-horizon discretization ----
    nfe_finite=2,
    ncp_finite=3,
    # ---- Infinite-horizon discretization ----
    nfe_infinite=3,
    ncp_infinite=3,
    # ---- Controller settings ----
    infinite_horizon=True,
    endpoint_constraints=True,
    custom_objective=True,
    input_suppression=True,
    input_suppression_factor=0.5e5,
    # ---- Cost function ----
    stage_cost_weights=[1e5, 1e3, 1e2, 1e2],
    beta=1.2,
    # ---- Output ----
    tee_flag=False,
    save_data=True,
    save_figure=True,
)

if __name__ == '__main__':
    mpc_loop(options)
