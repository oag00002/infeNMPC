"""
NMPC simulation options.

All fields have default values so that adding a new option in the future
never breaks existing code — callers that don't set the new field simply
get the default.
"""
from dataclasses import dataclass, field, replace
from typing import Any, List, Optional


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
    terminal_constraint_type : str
        Controls how the terminal state constraint is enforced.

        * ``'hard'`` — add equality constraints pinning each selected variable
          to its steady-state value at the terminal time.  For the infinite-
          horizon block this uses Lagrange extrapolation to τ=1 for algebraic
          variables (CVs without a ``DerivativeVar``, and all MVs), and direct
          access at τ=1 for differential state CVs.  For the finite-horizon
          block the last RADAU collocation point (= the FE right endpoint) is
          used directly.
        * ``'soft'`` — add a quadratic penalty term to the objective.  The
          per-variable weights come from ``stage_cost_weights``; an additional
          scalar ``terminal_soft_weight`` scales the total penalty.
        * ``'none'`` — no terminal constraint.

        Default: ``'hard'``.
    terminal_constraint_variables : str
        Which variables are included in the terminal constraint.

        * ``'cv'``   — ``CV_index`` only.
        * ``'mv'``   — ``MV_index`` only.
        * ``'cvmv'`` — both (default).
    terminal_soft_weight : float
        Scalar multiplier applied to the terminal soft-constraint penalty.
        Per-variable weights still come from ``stage_cost_weights``.
        Only used when ``terminal_constraint_type='soft'``.  Default: ``1.0``.
    objective : str
        Controls the MPC stage cost and steady-state operating point:

        * ``'economic'`` — use the model's ``custom_objective()`` function as
          the stage cost.  The steady-state operating point is found by
          minimising that economic cost.  Requires ``custom_objective`` to be
          defined in the model file.
        * ``'tracking'`` — use a quadratic tracking cost.  The setpoints used
          as tracking targets are controlled by ``tracking_setpoint``.

        Default: ``'economic'``.
    tracking_setpoint : str
        How setpoints for the quadratic tracking cost are determined.  Only
        relevant when ``objective='tracking'``.

        * ``'model'`` — use ``m.setpoints`` as declared in the model file.
          The steady-state operating point is found by minimising the distance
          to those setpoints.
        * ``'economic'`` — find the steady-state operating point by minimising
          the model's ``custom_objective()`` function, then use the resulting
          values as the tracking targets.  Requires ``custom_objective`` to be
          defined in the model file.

        Default: ``'model'``.
    terminal_cost_riemann : bool
        If True, approximate the terminal cost with a Riemann sum.
    input_suppression : bool
        If True, add a move-suppression penalty on MV increments.
    input_suppression_factor : float
        Weighting factor for the move-suppression penalty.
    stage_cost_weights : list of float
        Per-variable weights for the quadratic tracking stage cost.
    slack_penalty_weight : float
        Weight on soft-constraint slack variables (``m.slack_index``) in the
        tracking objective.  Added to the stage cost as
        ``slack_penalty_weight * Σ eps_var`` summed over all entries of every
        variable named in ``m.slack_index``.  Only applied when
        ``objective='tracking'``; ignored for ``objective='economic'`` since
        ``custom_objective`` is expected to include its own slack penalty.
        Set to 0 to disable.  Default: ``1e4``.
    gamma : float or None
        Time-compression parameter for the infinite-horizon transformation.
        If ``None`` (the default), gamma is chosen automatically so that the
        first collocation point in the infinite block satisfies
        ``tau_1 = tanh(gamma * sampling_time)``; i.e.
        ``gamma = atanh(tau_1) / sampling_time``.  Once computed it is stored
        on the infinite-horizon block as ``infinite_block.gamma``.
    beta : float
        Weighting factor on the terminal cost relative to the stage cost.
    initialization_assist : bool
        If True, warm-start via progressive sampling-time reduction.
    initialization_assist_sampling_time_start : float
        Starting (large) sampling time for the warm-start sequence.
    safe_run : bool
        If True, overwrite ``io_data.csv`` after every MPC iteration so that
        CV/MV data are preserved on solver failure or crash.  Slower than the
        default end-of-run save.  Use with the ``live_plot.py`` watcher for
        in-process monitoring.
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
    debug_flag : bool
        If True, print the most violated constraints after every controller
        solve (and before re-raising on solver failure).
    revive_run : int
        Maximum number of times to retry a failed controller solve using
        relaxed IPOPT settings (``bound_relax_factor=1e-8``, ``tol=1e-6``)
        before aborting the run.  ``0`` disables retry.  Not applied when the
        termination condition is ``infeasible`` or ``maxIterations``.
        Default: ``0``.
    dynamic_initial_conditions : bool
        If True, use the dynamic IC path: state variables are fixed to their
        ``initialize=`` values (with optional scalar overrides from
        ``initial_values``), and IPOPT solves only for algebraic variables and
        free derivatives. Derivatives are NOT fixed to zero, allowing the
        starting point to be off steady-state. If False (default), the existing
        steady-state NLP toward ``initial_values`` is used.
    model_output_dir : str or None
        If set, pprint each intermediate Pyomo model to a ``.txt`` file in
        this directory.  The following files are written:

        * ``ss_model.txt``      — steady-state setpoint model
        * ``iv_model.txt``      — initial-condition model
        * ``controller_model.txt`` — full controller NLP
        * ``plant_model.txt``   — plant integration model
        * ``dynamic_ic_model.txt`` — dynamic IC feasibility model (only when
          ``dynamic_initial_conditions=True``)

        The directory is created automatically if it does not exist.
        If ``None`` (default), no pprint files are written.
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
    terminal_constraint_type: str = 'hard'
    terminal_constraint_variables: str = 'cvmv'
    terminal_soft_weight: float = 1.0
    objective: str = 'economic'
    tracking_setpoint: str = 'model'
    initialize_with_initial_data: bool = False
    terminal_cost_riemann: bool = False
    initialization_assist: bool = False
    initialization_assist_sampling_time_start: float = 10.0

    # Input suppression
    input_suppression: bool = False
    input_suppression_factor: float = 1.0

    # Cost function parameters
    stage_cost_weights: List[float] = field(default_factory=lambda: [1.0, 1.0, 1 / 600])
    slack_penalty_weight: float = 1e4
    gamma: Optional[float] = None
    beta: float = 1.0

    # Lyapunov stability constraint
    lyap_flag: bool = False
    lyap_delta: float = 0.01
    lyap_constraint_type: str = 'hard'
    lyap_soft_weight: float = 1e4
    lyap_beta: float = 1.2

    # Display / data-output options
    safe_run: bool = False
    save_data: bool = True
    save_figure: bool = True

    # Disturbance options
    disturb_flag: bool = False
    disturb_distribution: str = 'normal'
    disturb_seeded: bool = True

    # Debugging
    debug_flag: bool = False
    revive_run: int = 0

    # Initial condition mode
    dynamic_initial_conditions: bool = False

    # Model inspection
    model_output_dir: Optional[str] = None

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
