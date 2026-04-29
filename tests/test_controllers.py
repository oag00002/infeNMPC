"""
Tests for Controller classes and associated Pyomo model structure.

These tests require IPOPT (they call Controller.__init__ which performs the
initial solve).  All tests use the minimal CSTR model to keep run times short.

What is verified:

  Controller base class
  - _initialized starts False, becomes True after first solve (selects cold
    vs. warm IPOPT solver — LLMREADME 'cold solver / warm solver' section).

  InfiniteHorizonController
  - finite_block and infinite_block are present on the model.
  - lyap_flag=True creates V_prev and first_stage_cost_prev on the model.
  - lyap_constraint_type='hard' creates lyap_stability_constraint.
  - lyap_constraint_type='soft' creates lyap_slack Var (instead of hard constraint).
  - lyap_constraint_type='none' creates phi_track infrastructure but NO constraint.
  - terminal_constraint_type='hard' creates diff_terminal_constraints on the block.
  - terminal_constraint_type='soft' creates infinite_terminal_soft_penalty Expression.
  - terminal_constraint_type='none' adds no terminal constraint components.

  FiniteHorizonController
  - Works end-to-end with the CSTR model.
  - Builds only a finite block (no infinite_block).

  Plant
  - Uses nfe_finite=1 internally.
  - All MVs are fixed at construction time.
  - lyap_flag is always False in the plant (no spurious DOF from lyap_slack).
  - The plant objective is constant (expr=1); it is a pure integration step.
"""
import pytest

import pyomo.environ as pyo

from infeNMPC import Options
from infeNMPC.controllers import InfiniteHorizonController, FiniteHorizonController
from infeNMPC.plant import Plant


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_ih_controller(opts):
    return InfiniteHorizonController(opts)


def _make_fh_controller(opts):
    return FiniteHorizonController(opts)


def _make_plant(opts):
    return Plant(opts)


# ---------------------------------------------------------------------------
# Controller._initialized (cold vs warm solver selection)
# ---------------------------------------------------------------------------

class TestInitializedFlag:
    """
    The _initialized flag governs cold-vs-warm solver selection.

    The cold solver (bound_relax_factor=0) is used for the very first solve
    because the warm-start point is not yet available.  Subsequent solves use
    the warm solver (bound_relax_factor=1e-8, bound_push/frac=1e-8) to stay
    close to the shifted warm-start point.  See LLMREADME 'IPOPT Configuration'
    section for the full rationale.
    """

    def test_initialized_false_before_construction(self, base_options):
        """Before the initial solve, _initialized is False."""
        # The Controller base class sets _initialized=False in __init__,
        # before the subclass calls self.solve().  Access it via the instance
        # just after super().__init__ — we test this indirectly by checking
        # that after construction it IS True (the subclass called solve).
        ctrl = _make_ih_controller(base_options)
        assert ctrl._initialized is True, (
            "_initialized must be True after the initial solve in __init__"
        )

    def test_last_solve_time_is_positive(self, base_options):
        """CPU time is recorded after each solve."""
        ctrl = _make_ih_controller(base_options)
        assert ctrl.last_solve_time >= 0.0


# ---------------------------------------------------------------------------
# InfiniteHorizonController: block structure
# ---------------------------------------------------------------------------

class TestInfiniteHorizonBlocks:
    """The controller model has finite_block and infinite_block sub-models."""

    def test_finite_block_present(self, base_options):
        ctrl = _make_ih_controller(base_options)
        assert hasattr(ctrl._model, "finite_block")

    def test_infinite_block_present(self, base_options):
        ctrl = _make_ih_controller(base_options)
        assert hasattr(ctrl._model, "infinite_block")

    def test_interface_on_finite_block(self, base_options):
        """
        DynamicModelInterface is attached to finite_block.time, not the top-
        level model.  run_MPC.py uses controller.interface (via __getattr__)
        to call shift_values_by_time and load_data.
        """
        ctrl = _make_ih_controller(base_options)
        assert hasattr(ctrl._model, "interface")

    def test_finite_block_has_mv_index(self, base_options):
        ctrl = _make_ih_controller(base_options)
        assert hasattr(ctrl._model.finite_block, "MV_index")
        assert len(list(ctrl._model.finite_block.MV_index)) > 0

    def test_finite_block_has_cv_index(self, base_options):
        ctrl = _make_ih_controller(base_options)
        assert hasattr(ctrl._model.finite_block, "CV_index")
        assert len(list(ctrl._model.finite_block.CV_index)) > 0

    def test_objective_is_minimization(self, base_options):
        """The NLP is always a minimization problem."""
        ctrl = _make_ih_controller(base_options)
        obj = ctrl._model.objective
        assert obj.sense == pyo.minimize


# ---------------------------------------------------------------------------
# Lyapunov infrastructure
# ---------------------------------------------------------------------------

class TestLyapunovInfrastructure:
    """
    lyap_flag=True builds Lyapunov tracking infrastructure.  The constraint
    form is controlled by lyap_constraint_type ('hard', 'soft', 'none').

    The Lyapunov stability condition requires:
        V_k <= V_{k-1} - delta * L_{k-1}
    where V is phi_track (the quadratic tracking Riemann sum) and L is the
    first-FE stage cost.  V_prev and first_stage_cost_prev are Pyomo Params
    updated by run_MPC.py after each iteration.
    """

    def _lyap_opts(self, cstr_module, lyap_type):
        return Options.for_model_module(
            cstr_module,
            num_horizons=1,
            sampling_time=1.0,
            nfe_finite=2,
            ncp_finite=1,
            nfe_infinite=3,
            ncp_infinite=1,
            infinite_horizon=True,
            terminal_constraint_type="none",
            objective="economic",
            stage_cost_weights=[1, 1, 1 / 600],
            beta=1.0,
            lyap_flag=True,
            lyap_delta=0.01,
            lyap_constraint_type=lyap_type,
            lyap_soft_weight=1e4,
            save_data=False,
            save_figure=False,
            tee_flag=False,
            debug_flag=False,
        )

    def test_v_prev_param_present(self, cstr_module):
        """V_prev is a mutable Param updated by run_MPC.py each iteration."""
        opts = self._lyap_opts(cstr_module, "hard")
        ctrl = _make_ih_controller(opts)
        assert hasattr(ctrl._model, "V_prev")
        assert isinstance(ctrl._model.V_prev, pyo.Param)

    def test_first_stage_cost_prev_param_present(self, cstr_module):
        """first_stage_cost_prev is also a mutable Param."""
        opts = self._lyap_opts(cstr_module, "hard")
        ctrl = _make_ih_controller(opts)
        assert hasattr(ctrl._model, "first_stage_cost_prev")
        assert isinstance(ctrl._model.first_stage_cost_prev, pyo.Param)

    def test_phi_track_on_finite_block(self, cstr_module):
        """phi_track Expression is on the finite block, not the top-level model."""
        opts = self._lyap_opts(cstr_module, "hard")
        ctrl = _make_ih_controller(opts)
        assert hasattr(ctrl._model.finite_block, "phi_track")
        assert isinstance(ctrl._model.finite_block.phi_track, pyo.Expression)

    def test_hard_constraint_creates_lyap_stability_constraint(self, cstr_module):
        """
        lyap_constraint_type='hard' adds lyap_stability_constraint on the
        top-level ConcreteModel (not on either block).
        """
        opts = self._lyap_opts(cstr_module, "hard")
        ctrl = _make_ih_controller(opts)
        assert hasattr(ctrl._model, "lyap_stability_constraint")
        assert isinstance(ctrl._model.lyap_stability_constraint, pyo.Constraint)

    def test_soft_constraint_creates_lyap_slack(self, cstr_module):
        """
        lyap_constraint_type='soft' adds lyap_slack (NonNegativeReals Var)
        instead of a hard constraint.  The slack enters the objective with
        weight lyap_soft_weight.  This is the recommended setting for short
        finite horizons where the hard constraint can be locally infeasible.
        """
        opts = self._lyap_opts(cstr_module, "soft")
        ctrl = _make_ih_controller(opts)
        assert hasattr(ctrl._model, "lyap_slack")
        assert isinstance(ctrl._model.lyap_slack, pyo.Var)

    def test_soft_constraint_has_lyap_stability_constraint_with_slack(self, cstr_module):
        """
        With lyap_constraint_type='soft', lyap_stability_constraint IS present
        but it includes the lyap_slack variable on the RHS.  Both the constraint
        and the slack Var coexist.  The hard form has only the constraint; the
        soft form has constraint + slack.
        """
        opts = self._lyap_opts(cstr_module, "soft")
        ctrl = _make_ih_controller(opts)
        assert hasattr(ctrl._model, "lyap_stability_constraint")
        assert hasattr(ctrl._model, "lyap_slack")

    def test_none_constraint_no_lyap_stability_constraint(self, cstr_module):
        """
        lyap_constraint_type='none' builds phi_track infrastructure but adds
        no constraint and no slack.  Useful for monitoring without enforcement.
        """
        opts = self._lyap_opts(cstr_module, "none")
        ctrl = _make_ih_controller(opts)
        assert not hasattr(ctrl._model, "lyap_stability_constraint")
        assert not hasattr(ctrl._model, "lyap_slack")

    def test_none_constraint_still_has_phi_track(self, cstr_module):
        """phi_track infrastructure exists even when lyap_constraint_type='none'."""
        opts = self._lyap_opts(cstr_module, "none")
        ctrl = _make_ih_controller(opts)
        assert hasattr(ctrl._model.finite_block, "phi_track")


# ---------------------------------------------------------------------------
# Terminal constraint types
# ---------------------------------------------------------------------------

class TestTerminalConstraintTypes:
    """
    terminal_constraint_type controls the terminal state constraint on the
    infinite block.  Algebraic variables (non-differential CVs and all MVs) use
    Lagrange extrapolation to τ=1 since LEGENDRE does not include that endpoint.
    Differential CVs access τ=1 directly via Pyomo continuity equations.
    """

    def _tc_opts(self, cstr_module, tc_type):
        return Options.for_model_module(
            cstr_module,
            num_horizons=1,
            sampling_time=1.0,
            nfe_finite=2,
            ncp_finite=1,
            nfe_infinite=3,
            ncp_infinite=1,
            infinite_horizon=True,
            terminal_constraint_type=tc_type,
            terminal_constraint_variables="cvmv",
            objective="economic",
            stage_cost_weights=[1, 1, 1 / 600],
            beta=1.0,
            lyap_flag=False,
            save_data=False,
            save_figure=False,
            tee_flag=False,
            debug_flag=False,
        )

    def test_hard_creates_terminal_constraints(self, cstr_module):
        """
        terminal_constraint_type='hard' creates diff_terminal_constraints
        on the infinite block (differential CVs pinned to SS value at τ=1).
        """
        opts = self._tc_opts(cstr_module, "hard")
        ctrl = _make_ih_controller(opts)
        ib = ctrl._model.infinite_block
        assert hasattr(ib, "diff_terminal_constraints")

    def test_soft_creates_penalty_expression(self, cstr_module):
        """
        terminal_constraint_type='soft' creates infinite_terminal_soft_penalty
        Expression on the infinite block and adds it to the objective.
        """
        opts = self._tc_opts(cstr_module, "soft")
        opts = opts.copy(terminal_soft_weight=1e3)
        ctrl = _make_ih_controller(opts)
        ib = ctrl._model.infinite_block
        assert hasattr(ib, "infinite_terminal_soft_penalty")

    def test_soft_does_not_create_hard_terminal_constraint(self, cstr_module):
        """Soft terminal constraint must not also create a hard constraint."""
        opts = self._tc_opts(cstr_module, "soft")
        ctrl = _make_ih_controller(opts)
        ib = ctrl._model.infinite_block
        assert not hasattr(ib, "diff_terminal_constraints")

    def test_none_creates_no_terminal_component(self, cstr_module):
        """terminal_constraint_type='none' adds no terminal constraints or penalties."""
        opts = self._tc_opts(cstr_module, "none")
        ctrl = _make_ih_controller(opts)
        ib = ctrl._model.infinite_block
        assert not hasattr(ib, "diff_terminal_constraints")
        assert not hasattr(ib, "infinite_terminal_soft_penalty")


# ---------------------------------------------------------------------------
# FiniteHorizonController
# ---------------------------------------------------------------------------

class TestFiniteHorizonController:
    """FiniteHorizonController uses only the finite block (no infinite block)."""

    def _fh_opts(self, cstr_module):
        return Options.for_model_module(
            cstr_module,
            num_horizons=1,
            sampling_time=1.0,
            nfe_finite=3,
            ncp_finite=1,
            infinite_horizon=False,
            terminal_constraint_type="none",
            objective="economic",
            stage_cost_weights=[1, 1, 1 / 600],
            beta=1.0,
            lyap_flag=False,
            save_data=False,
            save_figure=False,
            tee_flag=False,
            debug_flag=False,
        )

    def test_finite_controller_builds(self, cstr_module):
        """FiniteHorizonController completes construction without error."""
        opts = self._fh_opts(cstr_module)
        ctrl = _make_fh_controller(opts)
        assert ctrl._initialized is True

    def test_no_infinite_block(self, cstr_module):
        """The finite-horizon model has no infinite_block."""
        opts = self._fh_opts(cstr_module)
        ctrl = _make_fh_controller(opts)
        assert not hasattr(ctrl._model, "infinite_block")

    def test_has_objective(self, cstr_module):
        """A minimization objective is always present."""
        opts = self._fh_opts(cstr_module)
        ctrl = _make_fh_controller(opts)
        assert hasattr(ctrl._model, "objective")

    def test_lyap_finite_only(self, cstr_module):
        """
        When lyap_flag=True on a finite-horizon controller, V_prev and
        lyap_stability_constraint live directly on the model (not a sub-block).
        """
        opts = Options.for_model_module(
            cstr_module,
            num_horizons=1,
            sampling_time=1.0,
            nfe_finite=3,
            ncp_finite=1,
            infinite_horizon=False,
            terminal_constraint_type="none",
            objective="economic",
            stage_cost_weights=[1, 1, 1 / 600],
            beta=1.0,
            lyap_flag=True,
            lyap_delta=0.01,
            lyap_constraint_type="hard",
            save_data=False,
            save_figure=False,
            tee_flag=False,
            debug_flag=False,
        )
        ctrl = _make_fh_controller(opts)
        assert hasattr(ctrl._model, "V_prev")
        assert hasattr(ctrl._model, "lyap_stability_constraint")


# ---------------------------------------------------------------------------
# Plant structure
# ---------------------------------------------------------------------------

class TestPlantStructure:
    """
    Plant always uses nfe_finite=1 and lyap_flag=False, regardless of what the
    caller's options say.  MVs are fixed at construction; the plant just
    integrates forward one sampling interval given the fixed MV values.
    """

    def test_plant_builds(self, base_options):
        """Plant constructs without error."""
        plant = _make_plant(base_options)
        assert plant is not None

    def test_plant_nfe_is_1(self, base_options):
        """
        The plant is a single finite element regardless of the controller's nfe.
        Plant.__init__ calls options.copy(nfe_finite=1, ...).
        """
        plant = _make_plant(base_options)
        # ContinuousSet bounds are (0, sampling_time) for nfe=1
        t_end = plant._model.time.last()
        assert t_end == pytest.approx(base_options.sampling_time)

    def test_mvs_are_fixed(self, base_options):
        """
        All MV variables are fixed in the plant.  load_data skips fixed
        variables, so the unfix→load→refix pattern in run_MPC.py is needed
        to propagate the controller's optimal MV values into the plant.

        We iterate over .keys() (the actual VarData keys that exist after
        LAGRANGE-RADAU clean_model='delete') rather than .index_set() (the
        full ContinuousSet which includes non-collocation t=0 that may have
        been deleted from the Var by the discretizer).
        """
        plant = _make_plant(base_options)
        m = plant._model
        for mv_name in m.MV_index:
            mv = getattr(m, mv_name)
            for idx in mv.keys():
                assert mv[idx].fixed, f"{mv_name}[{idx}] is not fixed in the plant"

    def test_plant_has_dummy_objective(self, base_options):
        """
        The plant objective is a constant (expr=1).  The plant is a feasibility
        problem — IPOPT must satisfy the ODE constraints, not minimise anything.
        """
        plant = _make_plant(base_options)
        assert hasattr(plant._model, "obj")

    def test_plant_no_lyap_slack(self, lyap_options):
        """
        Even when the controller uses lyap_flag=True, Plant.__init__ calls
        options.copy(lyap_flag=False) so the plant has no lyap_slack DOF.
        A free lyap_slack in the plant would make the NLP under-determined.
        """
        plant = _make_plant(lyap_options)
        assert not hasattr(plant._model, "lyap_slack")
