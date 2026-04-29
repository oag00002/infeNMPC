"""
Tests for low-level model construction in infeNMPC/make_model.py.

These tests call _finite_block_gen and _infinite_block_gen directly on fresh
Pyomo ConcreteModel objects.  No IPOPT solve is required — the functions only
perform discretization and algebraic model setup.

Key invariants verified:

  1. t=0 constraint deletion (finite block)
     After equations_write runs, all constraint entries at t=0 are deleted.
     This is critical because LAGRANGE-RADAU does not use t=0 as a collocation
     point.  Constraints at t=0 couple algebraic CVs to state variables at a
     non-collocation time, causing incorrect ICs after iteration 0 when state
     vars are updated via shift_values_by_time.  The deletion was added in the
     April-2026 "Plant Squareness" fix session.

  2. Endpoint deletion in the infinite block (t=0 and t=1)
     LAGRANGE-LEGENDRE places collocation points strictly inside (0, 1).  Both
     endpoints are non-collocation points, so equations_write entries at both
     t=0 and t=1 must be deleted.  This mirrors the RADAU fix and was added
     as part of the same session.

  3. phi_track Expression when lyap_flag=True
     The Lyapunov Riemann-sum measure phi_track is built as a scalar Pyomo
     Expression on the finite block.  Verifying it exists prevents regressions
     to the old ODE-integral form which was replaced by the scalar Expression.

  4. state_vars and deriv_vars are stored on the block
     The framework relies on block.state_vars / block.deriv_vars throughout
     run_MPC.py for IC loading and Lyapunov updates.
"""
import sys
from pathlib import Path

import pyomo.environ as pyo
import pytest

# ---------------------------------------------------------------------------
# Make CSTR example importable
# ---------------------------------------------------------------------------
_CSTR_DIR = (
    Path(__file__).parent.parent
    / "examples"
    / "lyap_flag_examples"
    / "enmpc_cstr"
)
if str(_CSTR_DIR) not in sys.path:
    sys.path.insert(0, str(_CSTR_DIR))

import small_cstr_model  # noqa: E402

from infeNMPC import Options
from infeNMPC.make_model import _finite_block_gen, _infinite_block_gen
from pyomo.contrib.mpc import ScalarData


# ---------------------------------------------------------------------------
# CSTR steady-state values (approximate; not from a solve)
# These are needed by _build_phi_terms inside _finite_block_gen when
# lyap_flag=True, and by terminal constraint logic when terminal_constraint_type
# is not 'none'.  In the real flow, _make_finite_horizon_model obtains these
# from _solve_steady_state_model (an IPOPT solve).  Here we use approximate
# values so the block can be built without running IPOPT.
# ---------------------------------------------------------------------------
_CSTR_SS_VALS = {"Ca": 0.5, "Cb": 0.5, "Fa0": 12.0}
_CSTR_SS_DATA = ScalarData({"Ca[*]": 0.5, "Cb[*]": 0.5, "Fa0[*]": 12.0})


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def finite_opts():
    """
    Minimal Options for building a finite block without IPOPT.

    Uses terminal_constraint_type='none' and infinite_horizon=True so that
    _finite_block_gen does not try to add terminal constraints (which would
    need steady_state_values).  We pre-populate steady_state_values anyway
    in _build_finite_block so lyap tests also work.
    """
    return Options.for_model_module(
        small_cstr_model,
        nfe_finite=2,
        ncp_finite=3,
        sampling_time=1.0,
        infinite_horizon=True,   # prevents terminal constraints on this block
        terminal_constraint_type="none",
        lyap_flag=False,
        objective="economic",
        stage_cost_weights=[1, 1, 1 / 600],
    )


@pytest.fixture
def finite_lyap_opts():
    """Options for building a finite block with Lyapunov infrastructure."""
    return Options.for_model_module(
        small_cstr_model,
        nfe_finite=2,
        ncp_finite=1,
        sampling_time=1.0,
        infinite_horizon=True,   # prevents terminal constraints on this block
        terminal_constraint_type="none",
        lyap_flag=True,
        lyap_delta=0.01,
        objective="economic",
        stage_cost_weights=[1, 1, 1 / 600],
    )


@pytest.fixture
def infinite_opts():
    """Options for building an infinite block."""
    return Options.for_model_module(
        small_cstr_model,
        nfe_infinite=3,
        ncp_infinite=3,
        sampling_time=1.0,
        infinite_horizon=True,
        terminal_constraint_type="none",
        lyap_flag=False,
        objective="economic",
        stage_cost_weights=[1, 1, 1 / 600],
    )


def _build_finite_block(opts):
    """
    Build a fresh ConcreteModel and call _finite_block_gen.

    Pre-populates steady_state_values and steady_state_data so that phi_track
    expressions (lyap_flag=True) can be built without an IPOPT solve.
    """
    m = pyo.ConcreteModel()
    m.steady_state_values = _CSTR_SS_VALS
    m.steady_state_data = _CSTR_SS_DATA
    return _finite_block_gen(m, opts)


def _build_infinite_block(opts):
    """Build a fresh ConcreteModel and call _infinite_block_gen."""
    m = pyo.ConcreteModel()
    m.steady_state_values = _CSTR_SS_VALS
    m.steady_state_data = _CSTR_SS_DATA
    m.ss_obj_value = pyo.Param(initialize=0.0)
    return _infinite_block_gen(m, opts)


# ---------------------------------------------------------------------------
# Finite block: t=0 constraint deletion
# ---------------------------------------------------------------------------

class TestFiniteBlockT0Deletion:
    """
    After _finite_block_gen runs, no constraint should have an entry at t=0.

    The deletion is required because:
      - equations_write runs post-discretization (only then is m.time populated).
      - This causes every Constraint(m.time, rule=...) to create an entry at t=0.
      - t=0 is NOT a RADAU collocation point.
      - If left in, entries at t=0 couple algebraic CVs to differential state vars
        at a non-integrated time point, corrupting ICs after the first MPC shift.
    """

    def test_no_constraints_at_t0(self, finite_opts):
        """Zero constraint entries at t=0 after finite block generation."""
        m = _build_finite_block(finite_opts)
        t0 = m.time.first()
        violations = []
        for con in m.component_objects(pyo.Constraint, active=True):
            for idx in con.keys():
                t = idx[-1] if isinstance(idx, tuple) else idx
                if t == t0:
                    violations.append((con.name, idx))
        assert violations == [], (
            f"Found {len(violations)} constraint(s) at t=0 that should have been deleted: "
            + str(violations[:5])
        )

    def test_constraints_exist_at_nonzero_times(self, finite_opts):
        """
        Verify the deletion is selective — constraints at non-zero collocation
        points must still be present.  If all constraints were deleted, the
        model would be degenerate and this sanity check catches that.
        """
        m = _build_finite_block(finite_opts)
        t0 = m.time.first()
        count = sum(
            1
            for con in m.component_objects(pyo.Constraint, active=True)
            for idx in con.keys()
            if (idx[-1] if isinstance(idx, tuple) else idx) != t0
        )
        assert count > 0, "All constraints were deleted, including non-t0 entries"

    def test_state_vars_stored_on_block(self, finite_opts):
        """block.state_vars is populated — used by run_MPC.py for IC loading."""
        m = _build_finite_block(finite_opts)
        assert hasattr(m, "state_vars")
        assert len(m.state_vars) > 0

    def test_deriv_vars_stored_on_block(self, finite_opts):
        """block.deriv_vars is populated — used by the framework to identify ODEs."""
        m = _build_finite_block(finite_opts)
        assert hasattr(m, "deriv_vars")
        assert len(m.deriv_vars) > 0

    def test_mv_index_present(self, finite_opts):
        """MV_index is a non-empty Set after variables_initialize."""
        m = _build_finite_block(finite_opts)
        assert hasattr(m, "MV_index")
        assert len(list(m.MV_index)) > 0

    def test_cv_index_present(self, finite_opts):
        """CV_index is a non-empty Set after variables_initialize."""
        m = _build_finite_block(finite_opts)
        assert hasattr(m, "CV_index")
        assert len(list(m.CV_index)) > 0


# ---------------------------------------------------------------------------
# Finite block: Lyapunov infrastructure
# ---------------------------------------------------------------------------

class TestFiniteBlockLyapunov:
    """
    When lyap_flag=True, _finite_block_gen adds a scalar phi_track Expression.

    phi_track is a Riemann sum of the quadratic tracking cost at each finite-
    element right-endpoint.  It replaced the old ODE-integral form (a Var with
    DerivativeVar and phi_track_ode Constraint).  The Expression form is cheaper
    to evaluate and has no ODEs to integrate.
    """

    def test_phi_track_expression_present(self, finite_lyap_opts):
        """phi_track is a pyo.Expression on the block."""
        m = _build_finite_block(finite_lyap_opts)
        assert hasattr(m, "phi_track")
        assert isinstance(m.phi_track, pyo.Expression)

    def test_phi_track_not_a_var(self, finite_lyap_opts):
        """phi_track must not be a Var (old ODE-integral form was a Var)."""
        m = _build_finite_block(finite_lyap_opts)
        assert not isinstance(m.phi_track, pyo.Var)

    def test_phi_track_evaluates_to_nonnegative(self, finite_lyap_opts):
        """
        phi_track is a sum of squared terms — it must always be >= 0.
        The initial value (default initialize=) may not be at the setpoint, so
        the value could be > 0, but never negative.
        """
        m = _build_finite_block(finite_lyap_opts)
        val = pyo.value(m.phi_track)
        assert val >= 0.0

    def test_phi_track_absent_when_lyap_disabled(self, finite_opts):
        """phi_track is not added when lyap_flag=False."""
        m = _build_finite_block(finite_opts)
        assert not hasattr(m, "phi_track")


# ---------------------------------------------------------------------------
# Infinite block: endpoint deletion (t=0 and t=1)
# ---------------------------------------------------------------------------

class TestInfiniteBlockEndpointDeletion:
    """
    LAGRANGE-LEGENDRE places collocation points strictly inside (0, 1).
    Both endpoints are non-collocation: equations_write entries at t=0 and t=1
    must be deleted.

    This fix mirrors the RADAU fix and was applied in the same session.
    The infinite-block deletion is critical for the same reason: spurious
    constraints at t=0/t=1 couple variables at non-integrated times, causing
    IC corruption and warm-start failures.
    """

    def test_no_constraints_at_t0(self, infinite_opts):
        """Zero constraint entries at τ=0 after infinite block generation."""
        m = _build_infinite_block(infinite_opts)
        t0 = m.time.first()
        violations = []
        for con in m.component_objects(pyo.Constraint, active=True):
            for idx in con.keys():
                t = idx[-1] if isinstance(idx, tuple) else idx
                if t == t0:
                    violations.append((con.name, idx))
        assert violations == [], (
            f"Found {len(violations)} constraint(s) at τ=0 in infinite block: "
            + str(violations[:5])
        )

    def test_no_constraints_at_t1(self, infinite_opts):
        """
        Zero constraint entries at τ=1 after infinite block generation.
        τ=1 is the terminal point of the transformed domain and is also not
        a LEGENDRE collocation point.  Constraints there would incorrectly
        link algebraic terminal variables to collocation-point physics.
        """
        m = _build_infinite_block(infinite_opts)
        t1 = m.time.last()
        violations = []
        for con in m.component_objects(pyo.Constraint, active=True):
            for idx in con.keys():
                t = idx[-1] if isinstance(idx, tuple) else idx
                if t == t1:
                    violations.append((con.name, idx))
        # terminal_cost and phi-ODE constraints are explicitly skipped in
        # _add_dphidt (they rule on t==0 and t==1 internally), so they should
        # not appear at the endpoints.
        assert violations == [], (
            f"Found {len(violations)} constraint(s) at τ=1 in infinite block: "
            + str(violations[:5])
        )

    def test_constraints_exist_at_interior_collocation_points(self, infinite_opts):
        """
        Sanity check: constraints at interior collocation points must survive.
        The deletion is selective — it must not remove all constraints.
        """
        m = _build_infinite_block(infinite_opts)
        t0, t1 = m.time.first(), m.time.last()
        count = sum(
            1
            for con in m.component_objects(pyo.Constraint, active=True)
            for idx in con.keys()
            if (idx[-1] if isinstance(idx, tuple) else idx) not in (t0, t1)
        )
        assert count > 0, "No constraints remain at interior collocation points"
