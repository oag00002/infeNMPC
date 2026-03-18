"""
NMPC simulation options.

All fields have default values so that adding a new option in the future
never breaks existing code — callers that don't set the new field simply
get the default.
"""
from dataclasses import dataclass, field, replace
from typing import Any, List


@dataclass
class Options:
    """
    Configuration for infinite and finite horizon NMPC simulations.

    Every attribute has a default value so that constructing ``Options()``
    is always valid and new fields can be added without breaking callers.

    Attributes
    ----------
    model_name : str
        Name of the model module in ``infeNMPC/models/``.
    num_horizons : int
        Number of MPC iterations to run.
    nfe_finite : int
        Number of finite elements in the finite-horizon block.
    ncp_finite : int
        Number of collocation points per finite-horizon finite element.
    sampling_time : float
        Time interval between successive control updates.
    infinite_horizon : bool
        If True, use infinite-horizon MPC; otherwise finite-horizon.
    nfe_infinite : int
        Number of finite elements for the infinite-horizon approximation.
    ncp_infinite : int
        Number of collocation points per infinite-horizon finite element.
    tee_flag : bool
        If True, print IPOPT output to the console.
    endpoint_constraints : bool
        If True, enforce a terminal state equality constraint at t=1 on the
        infinite-horizon block.
    custom_objective : bool
        If True, use the model's ``custom_objective`` economic stage cost.
    terminal_cost_riemann : bool
        If True, approximate the terminal cost with a Riemann sum.
    input_suppression : bool
        If True, add a move-suppression penalty on MV increments.
    input_suppression_factor : float
        Weighting factor for the move-suppression penalty.
    stage_cost_weights : list of float
        Per-variable weights for the quadratic tracking stage cost.
    gamma : float
        Time-compression parameter for the infinite-horizon transformation.
    beta : float
        Weighting factor on the terminal cost relative to the stage cost.
    initialization_assist : bool
        If True, warm-start via progressive sampling-time reduction.
    initialization_assist_sampling_time_start : float
        Starting (large) sampling time for the warm-start sequence.
    live_plot : bool
        If True, update a live plot at every MPC iteration.
    plot_end : bool
        If True, display a summary plot after the run.
    save_data : bool
        If True, write result CSVs to disk.
    save_figure : bool
        If True, save the summary figure to disk.
    disturb_flag : bool
        If True, apply random disturbances to the plant.
    disturb_distribution : str
        Distribution type for disturbances (e.g. ``'normal'``).
    disturb_seeded : bool
        If True, use a fixed random seed for reproducibility.
    """

    # Model selection
    # Supply *either* model_module (a live Python module object) *or* model_name
    # (a dotted-importable name).  model_module takes priority when set.
    model_module: Any = field(default=None, repr=False, compare=False)
    model_name: str = 'unknown'

    # Simulation control
    num_horizons: int = 100
    nfe_finite: int = 2
    ncp_finite: int = 3
    sampling_time: float = 1.0

    # Infinite horizon settings
    infinite_horizon: bool = True
    nfe_infinite: int = 3
    ncp_infinite: int = 3

    # Solver and model options
    tee_flag: bool = False
    endpoint_constraints: bool = True
    custom_objective: bool = True
    initialize_with_initial_data: bool = False
    terminal_cost_riemann: bool = False
    initialization_assist: bool = False
    initialization_assist_sampling_time_start: float = 10.0

    # Input suppression
    input_suppression: bool = False
    input_suppression_factor: float = 1.0

    # Cost function parameters
    stage_cost_weights: List[float] = field(default_factory=lambda: [1.0, 1.0, 1 / 600])
    gamma: float = 0.0375847
    beta: float = 1.0

    # Display / data-output options
    live_plot: bool = False
    plot_end: bool = True
    save_data: bool = True
    save_figure: bool = True

    # Disturbance options
    disturb_flag: bool = False
    disturb_distribution: str = 'normal'
    disturb_seeded: bool = True

    def copy(self, **overrides) -> 'Options':
        """Return a shallow copy of this Options with selected fields overridden."""
        return replace(self, **overrides)

    @classmethod
    def for_model_module(cls, module, **overrides) -> 'Options':
        """
        Create an Options pre-configured from a live module object.

        Applies the module's ``default_options()`` dict (if present) on top of
        the global defaults, then applies any explicit *overrides*.  Sets
        ``model_module`` to *module* and derives ``model_name`` from the
        module's ``__name__`` attribute.

        Parameters
        ----------
        module
            A Python module implementing ``variables_initialize``,
            ``equations_write``, and optionally ``custom_objective`` and
            ``default_options``.
        **overrides
            Any ``Options`` fields to override after applying model defaults.
        """
        raw_name = getattr(module, '__name__', None) or 'custom'
        name = raw_name.split('.')[-1]
        instance = cls(model_module=module, model_name=name)
        if hasattr(module, 'default_options'):
            for k, v in module.default_options().items():
                setattr(instance, k, v)
        for k, v in overrides.items():
            setattr(instance, k, v)
        return instance

    @classmethod
    def for_model(cls, model_name: str, **overrides) -> 'Options':
        """
        Create an Options pre-configured for a model installed as a Python package.

        Applies the model's ``default_options()`` dict (if present) on top of
        the global defaults, then applies any explicit *overrides*.  The model
        is loaded via ``_load_model`` using ``model_name`` as a package-relative
        name (e.g. a module in ``infeNMPC.models``).

        For models that live in a local file, prefer ``for_model_module``
        instead.
        """
        from .model_equations import _load_model
        mod = _load_model(model_name)
        return cls.for_model_module(mod, model_name=model_name, **overrides)


def _import_settings() -> Options:
    """
    Return a default ``Options`` instance.

    Kept for backwards compatibility — new code should construct
    ``Options()`` or ``Options.for_model(name)`` directly.
    """
    return Options()
