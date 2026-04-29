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
    nfe_finite=10,
    ncp_finite=3,
    # ---- Infinite-horizon discretization ----
    infinite_horizon=True,
    nfe_infinite=3,
    ncp_infinite=3,
    # ---- Controller settings ----
    terminal_constraint_type='soft',
    terminal_soft_weight=1e4,
    objective='economic',
    # ---- Cost function ----
    stage_cost_weights=[1e5, 1e3, 1e2, 1e2],
    beta=1.0,
    # ---- Output ----
    tee_flag=True,
    safe_run=True,
    save_data=True,
    save_figure=True,
    # ---- Lyapunov stability constraint ----
    lyap_flag=True,
    lyap_delta=0.5,
    lyap_constraint_type='soft',
    lyap_soft_weight=1e4,   # rho in AMPL: penalises desceps
    lyap_beta=1.2,
    # ---- Debugging ----
    debug_flag=True,
    revive_run=1,
)

if __name__ == '__main__':
    mpc_loop(options)
