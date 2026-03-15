"""
infeNMPC - Infinite Horizon Nonlinear Model Predictive Control

A framework for running infinite and finite horizon NMPC simulations
on dynamic systems using Pyomo/IDAES. The core workflow is:

  1. Select (or write) a model module in ``infeNMPC/models/``.
  2. Configure simulation parameters via the ``Options`` class.
  3. Run the closed-loop MPC simulation with ``mpc_loop``.

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
infNMPC_options
    Simulation and solver configuration (the ``Options`` dataclass).
model_equations
    Model loading utilities (``_get_model``, ``_load_model``).
indexing_tools
    Utilities for variable indexing and Pyomo expression construction.
initialization_tools
    Progressive sampling-time reduction for robust controller warm-starting.
data_save_and_plot
    Result serialization (CSV), live plotting, and post-run figure generation.

Typical usage
-------------
Write a ``model.py`` implementing ``variables_initialize``, ``equations_write``,
and optionally ``custom_objective``.  Then in a ``run.py``:

    import model
    from infeNMPC import Options, mpc_loop

    options = Options.for_model_module(
        model,
        num_horizons=100,
        sampling_time=1.0,
        ...
    )
    mpc_loop(options)

See the ``examples/`` directory for complete worked examples.
"""

from .infNMPC_options import Options, _import_settings
from .plant import Plant
from .controllers import InfiniteHorizonController, FiniteHorizonController
from .run_MPC import mpc_loop

__all__ = [
    "Options",
    "_import_settings",
    "Plant",
    "InfiniteHorizonController",
    "FiniteHorizonController",
    "mpc_loop",
]
