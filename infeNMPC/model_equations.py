"""
Model loading utilities.

Two mechanisms are supported:

1. **Module object** (primary): set ``options.model_module`` to a Python module
   that implements ``variables_initialize``, ``equations_write``, and
   optionally ``custom_objective`` and ``default_options``.  This is how the
   examples in ``examples/`` work — each ``run.py`` imports its local
   ``model.py`` and passes it via ``Options.for_model_module(model)``.

2. **Package import by name** (legacy): set ``options.model_name`` to a
   module name importable from ``infeNMPC.models`` (or any other installed
   package).  ``_load_model`` handles this path.

All framework internals call ``_get_model(options)`` which automatically
selects the right mechanism.
"""
from importlib import import_module as _import_module


def _load_model(model_name: str):
    """
    Import and return a model module by name.

    The module must be importable as ``infeNMPC.models.<model_name>`` or as a
    top-level package named ``model_name``.

    Parameters
    ----------
    model_name : str
        Module name relative to ``infeNMPC.models``, e.g. ``'enmpc_cstr'``.

    Returns
    -------
    module
    """
    _MODELS_PACKAGE = 'infeNMPC.models'
    return _import_module(f'.{model_name}', package=_MODELS_PACKAGE)


def _get_model(options):
    """
    Return the model module for *options*.

    Uses ``options.model_module`` directly when set; otherwise falls back to
    loading by ``options.model_name`` via ``_load_model``.

    Parameters
    ----------
    options : Options

    Returns
    -------
    module
    """
    if options.model_module is not None:
        return options.model_module
    return _load_model(options.model_name)
