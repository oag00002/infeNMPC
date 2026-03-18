"""
Run script for the non-isothermal CSTR example.

System:  Non-isothermal CSTR — exothermic reaction A+B → C
States:  Ca, Cb, Cc, Cm (mol/L), T (K)
MVs:     Fa0 — feed rate of A (mol/min), mc — coolant flow rate (mol/min)
Objective: quadratic tracking of Cc and T to setpoints (5.18, 396.1)

Usage
-----
From the repository root:

    python examples/nonisothermal_cstr/run.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

import model  # noqa: E402

from infeNMPC import Options, mpc_loop  # noqa: E402

options = Options.for_model_module(
    model,
    # ---- Simulation control ----
    num_horizons=400,
    sampling_time=0.0025,
    # ---- Finite-horizon discretization ----
    nfe_finite=50,
    ncp_finite=5,
    # ---- Infinite-horizon discretization ----
    nfe_infinite=10,
    ncp_infinite=5,
    # ---- Controller settings ----
    infinite_horizon=True,
    endpoint_constraints=True,
    custom_objective=False,
    initialize_with_initial_data=True,
    # ---- Cost function ----
    stage_cost_weights=[1, 1e-2, 1e-2, 1e-3],
    gamma=0.015,
    beta=1.2,
    # ---- Output ----
    tee_flag=False,
    save_data=True,
    save_figure=True,
    plot_end=True,
)

if __name__ == '__main__':
    mpc_loop(options)
