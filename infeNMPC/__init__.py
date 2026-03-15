"""
infeNMPC - Infinite Horizon Nonlinear Model Predictive Control

A framework for running infinite and finite horizon NMPC simulations
on dynamic systems using Pyomo/IDAES. The core workflow is:

  1. Select (or write) a model module in ``infeNMPC/models/``.
  2. Configure simulation parameters via the ``Options`` class.
  3. Run the closed-loop MPC simulation with ``_mpc_loop``.

Modules
-------
run_MPC
    Main MPC simulation loop entry point.
controllers
    ``InfiniteHorizonController`` and ``FiniteHorizonController`` classes.
plant
    ``Plant`` class for the closed-loop dynamic simulation model.
make_model
    Pyomo model construction: steady-state, finite-horizon, infinite-horizon.
models/
    Per-system model modules (variables, equations, optional custom objective).
infNMPC_options
    Simulation and solver configuration (the ``Options`` dataclass).
indexing_tools
    Utilities for variable indexing and Pyomo expression construction.
initialization_tools
    Progressive sampling-time reduction for robust controller warm-starting.
data_save_and_plot
    Result serialization (CSV), live plotting, and post-run figure generation.
"""

from .infNMPC_options import Options, _import_settings
from .plant import Plant
from .controllers import InfiniteHorizonController, FiniteHorizonController
from .run_MPC import _mpc_loop

__all__ = [
    "Options",
    "_import_settings",
    "Plant",
    "InfiniteHorizonController",
    "FiniteHorizonController",
    "_mpc_loop",
]
