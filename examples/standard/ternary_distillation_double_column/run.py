"""
Run script for the eNMPC ternary distillation two-column example.

System:   Two distillation columns in series separating A, B, C
States:   x1[tray,comp], M1[tray] (Col 1); x2[tray,comp], M2[tray] (Col 2) — 246 total
MVs:      VB1, LT1, D1, B1 (Column 1); VB2, LT2, D2, B2 (Column 2)
CVs:      xD1A — A purity in Col 1 distillate
          xD2B — B purity in Col 2 distillate
          xC   — C purity in Col 2 bottoms
Objective: minimize feed + energy cost minus product revenue (economic NMPC)

Usage
-----
From the repository root:

    python examples/standard/ternary_distillation_double_column/run.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

import model  # noqa: E402

from infeNMPC import Options, mpc_loop  # noqa: E402

options = Options.for_model_module(
    model,
    # ---- Simulation control ----
    num_horizons=10,
    sampling_time=1,
    # ---- Finite-horizon discretization ----
    nfe_finite=2,
    ncp_finite=3,
    # ---- Infinite-horizon discretization ----
    infinite_horizon=True,
    nfe_infinite=3,
    ncp_infinite=3,
    # ---- Controller settings ----
    endpoint_constraints=False,
    custom_objective=True,
    terminal_cost_riemann=False,
    initialize_with_initial_data=False,
    initialization_assist=False,
    input_suppression=True,
    input_suppression_factor=1.0e3,
    # ---- Cost function ----
    # Weights order: [xD1A, xD2B, xC, VB1, LT1, D1, B1, VB2, LT2, D2, B2]
    stage_cost_weights=[1.0e4, 1.0e4, 1.0e4, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    gamma=0.05,
    beta=1.0,
    # ---- Output ----
    tee_flag=True,
    save_data=True,
    save_figure=True,
    # ---- Disturbances ----
    disturb_flag=False,
    disturb_distribution='normal',
    disturb_seeded=True,
)

if __name__ == '__main__':
    mpc_loop(options)
