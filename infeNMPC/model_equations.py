"""
Model equation dispatcher.

Loads a model module from ``infeNMPC/models/`` by name and provides
``_load_model`` for use by the rest of the framework.  Each model module
must implement:

    variables_initialize(m) -> m
    equations_write(m) -> m
    custom_objective(m, options) -> callable   [if options.custom_objective]
    default_options() -> dict                  [optional convenience]
"""
from importlib import import_module as _import_module

_MODELS_PACKAGE = 'infeNMPC.models'


def _load_model(model_name: str):
    """
    Import and return the model module identified by ``model_name``.

    Parameters
    ----------
    model_name : str
        The short name of the model (e.g. ``'enmpc_cstr'``, ``'pendulum'``).
        Must match a module inside the ``infeNMPC/models/`` package.

    Returns
    -------
    module
        The imported model module exposing ``variables_initialize``,
        ``equations_write``, and optionally ``custom_objective``.

    Raises
    ------
    ModuleNotFoundError
        If no module named ``model_name`` exists in ``infeNMPC/models/``.
    """
    return _import_module(f'.{model_name}', package=_MODELS_PACKAGE)
