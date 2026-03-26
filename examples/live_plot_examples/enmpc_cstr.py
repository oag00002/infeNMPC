"""
eNMPC CSTR example with live plotting in a single script.

Starts the live plotter as a background process, runs the MPC loop with
``safe_run=True`` so that ``io_data.csv`` is updated after every iteration,
then shuts the plotter down when the run finishes.

Usage
-----
From the repository root::

    python examples/live_plot_examples/enmpc_cstr.py
"""
import sys
from pathlib import Path

# Resolve the enmpc_cstr model relative to this file.
sys.path.insert(0, str(Path(__file__).parent.parent / 'standard' / 'enmpc_cstr'))
import model  # noqa: E402  (local model.py)

from infeNMPC import Options, mpc_loop                          # noqa: E402
from infeNMPC.live_plot import start_live_plotter, show_final_results  # noqa: E402

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
    safe_run=True,   # write io_data.csv after every iteration
)

if __name__ == '__main__':
    plotter = start_live_plotter(options)
    try:
        mpc_loop(options)
    finally:
        plotter.terminate()
        plotter.wait()
    show_final_results(options)
