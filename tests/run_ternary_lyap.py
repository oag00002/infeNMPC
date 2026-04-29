"""
Manual inspection test for the lyap-flag ternary distillation example.

Runs a short (2-iteration) closed-loop simulation and saves pprint snapshots
of every intermediate Pyomo model to tests/outputs/ so they can be inspected
by hand.  The output files are untracked (see .gitignore in this directory).

Files written to tests/outputs/
--------------------------------
ss_model.txt         -- steady-state setpoint NLP (finds operating point)
iv_model.txt         -- initial-condition NLP    (finds consistent t=0 state)
controller_model.txt -- full infinite-horizon controller NLP
plant_model.txt      -- single-step plant integration model

Usage (from repo root)
----------------------
    python tests/run_ternary_lyap.py
"""
import sys
from pathlib import Path

# Allow importing the example model from its own directory
EXAMPLE_DIR = Path(__file__).parent.parent / "examples" / "lyap_flag_examples" / "enmpc_ternary_distillation"
sys.path.insert(0, str(EXAMPLE_DIR))

import ternary_distillation_model  # noqa: E402

from infeNMPC import Options, mpc_loop  # noqa: E402

OUTPUT_DIR = str(Path(__file__).parent / "outputs")

options = Options.for_model_module(
    ternary_distillation_model,
    # ---- Simulation control (short run for inspection) ----
    num_horizons=2,
    sampling_time=1,
    # ---- Finite-horizon discretization ----
    nfe_finite=3,
    ncp_finite=3,
    # ---- Infinite-horizon discretization ----
    infinite_horizon=True,
    nfe_infinite=3,
    ncp_infinite=3,
    # ---- Controller settings ----
    terminal_constraint_type='hard',
    objective='economic',
    terminal_cost_riemann=False,
    initialize_with_initial_data=False,
    initialization_assist=False,
    input_suppression=True,
    input_suppression_factor=1.0e3,
    # ---- Cost function ----
    # Weights order: [xD1A, xD2B, xC, VB1, LT1, D1, B1, VB2, LT2, D2, B2]
    stage_cost_weights=[1.0e4, 1.0e4, 1.0e4, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    beta=1.0,
    # ---- Output ----
    tee_flag=False,
    safe_run=False,
    save_data=False,
    save_figure=False,
    # ---- Disturbances ----
    disturb_flag=False,
    disturb_distribution='normal',
    disturb_seeded=True,
    # ---- Lyapunov stability constraint ----
    lyap_flag=True,
    lyap_delta=0.01,
    # ---- Debugging ----
    debug_flag=True,
    # ---- Model inspection outputs ----
    model_output_dir=OUTPUT_DIR,
)

if __name__ == '__main__':
    mpc_loop(options)
