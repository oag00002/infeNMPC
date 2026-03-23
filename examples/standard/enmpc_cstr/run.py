"""
Run script for the eNMPC CSTR example.

System:  Isothermal CSTR — consecutive reaction A → B
States:  Ca (mol/L), Cb (mol/L)
MV:      Fa0 — feed rate of A (kmol/h)
Objective: maximise B production (economic infinite-horizon)

Usage
-----
From the repository root:

    python examples/enmpc_cstr/run.py

Or via the infeNMPC CLI (uses model defaults, no run-file customisation):

    python -m infeNMPC --model examples/enmpc_cstr/model.py
"""
import sys
from pathlib import Path

# Allow plain `import model` to resolve to model.py in this directory.
sys.path.insert(0, str(Path(__file__).parent))

import model  # noqa: E402  (local model.py)

from infeNMPC import Options, mpc_loop  # noqa: E402

options = Options.for_model_module(
    model,
    # ---- Simulation control ----
    num_horizons=100,
    sampling_time=1,
    # ---- Finite-horizon discretization ----
    nfe_finite=2,
    ncp_finite=1,
    # ---- Infinite-horizon discretization ----
    nfe_infinite=5,
    ncp_infinite=1,
    # ---- Controller settings ----
    infinite_horizon=True,
    endpoint_constraints=False,
    custom_objective=True,
    # ---- Cost function ----
    stage_cost_weights=[1, 1, 1 / 600],
    gamma=0.001,
    beta=1,
    # ---- Output ----
    tee_flag=False,
    save_data=True,
    save_figure=True,
)

if __name__ == '__main__':
    mpc_loop(options)
