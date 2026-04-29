"""
Tests for infeNMPC/tools/solver_tools.py.

solver_tools.py contains two utilities added in the April-2026 solver-hardening
session:

  _clip_to_bounds(model)
      Projects all variable values to [lb, ub] after each IPOPT solve.
      This is necessary because bound_relax_factor=1e-8 in the warm solver
      allows IPOPT to return values that are slightly outside declared bounds.
      Without clipping, downstream Pyomo set_value() calls on NonNegativeReals
      slacks would raise ValueError when the value is a small negative number.

  _attempt_revive(controller, i, options, caught_exc)
      Retries a failed controller solve up to options.revive_run times using
      a more tolerant IPOPT configuration.  Skips retry for infeasible /
      maxIterations termination conditions (where retrying is pointless).

These tests use simple hand-built Pyomo models to avoid requiring IPOPT.
"""
import pytest
import pyomo.environ as pyo

from infeNMPC.tools.solver_tools import _clip_to_bounds


# ---------------------------------------------------------------------------
# _clip_to_bounds
# ---------------------------------------------------------------------------

class TestClipToBounds:
    """_clip_to_bounds clips variable values to [lb, ub]."""

    def _make_model(self, val, lb=None, ub=None):
        m = pyo.ConcreteModel()
        m.x = pyo.Var(bounds=(lb, ub), initialize=val)
        return m

    def test_value_below_lb_is_clipped(self):
        """A value just below lb (e.g., -1e-9 for NonNegativeReals) is set to lb."""
        m = self._make_model(val=-1e-9, lb=0.0, ub=None)
        _clip_to_bounds(m)
        assert pyo.value(m.x) == 0.0

    def test_value_above_ub_is_clipped(self):
        """A value just above ub is set to ub."""
        m = self._make_model(val=1.0 + 1e-9, lb=0.0, ub=1.0)
        _clip_to_bounds(m)
        assert pyo.value(m.x) == pytest.approx(1.0)

    def test_value_in_range_unchanged(self):
        """A value already within [lb, ub] must not be modified."""
        m = self._make_model(val=0.5, lb=0.0, ub=1.0)
        _clip_to_bounds(m)
        assert pyo.value(m.x) == pytest.approx(0.5)

    def test_no_lower_bound(self):
        """Variables with lb=None are not clipped from below."""
        m = self._make_model(val=-100.0, lb=None, ub=0.0)
        _clip_to_bounds(m)
        assert pyo.value(m.x) == pytest.approx(-100.0)

    def test_no_upper_bound(self):
        """Variables with ub=None are not clipped from above."""
        m = self._make_model(val=100.0, lb=0.0, ub=None)
        _clip_to_bounds(m)
        assert pyo.value(m.x) == pytest.approx(100.0)

    def test_no_bounds_unchanged(self):
        """Variables with no bounds at all are not touched."""
        m = self._make_model(val=-999.0, lb=None, ub=None)
        _clip_to_bounds(m)
        assert pyo.value(m.x) == pytest.approx(-999.0)

    def test_handles_none_value(self):
        """Variables with value=None (uninitialised) are skipped without error."""
        m = pyo.ConcreteModel()
        m.x = pyo.Var(bounds=(0, 1))  # no initialize → value is None
        _clip_to_bounds(m)  # must not raise

    def test_indexed_var_all_indices_clipped(self):
        """_clip_to_bounds processes every index of an indexed Var."""
        m = pyo.ConcreteModel()
        m.idx = pyo.Set(initialize=[1, 2, 3])
        m.x = pyo.Var(m.idx, bounds=(0.0, 1.0))
        m.x[1].set_value(-0.001)   # below lb
        m.x[2].set_value(0.5)      # in range
        m.x[3].set_value(1.001)    # above ub
        _clip_to_bounds(m)
        assert pyo.value(m.x[1]) == pytest.approx(0.0)
        assert pyo.value(m.x[2]) == pytest.approx(0.5)
        assert pyo.value(m.x[3]) == pytest.approx(1.0)

    def test_nested_block_clipped(self):
        """
        The controller model uses nested Blocks (finite_block, infinite_block).
        _clip_to_bounds must descend into sub-blocks via descend_into=True.
        """
        m = pyo.ConcreteModel()
        m.sub = pyo.Block()
        m.sub.y = pyo.Var(bounds=(0, None), initialize=-0.5)
        _clip_to_bounds(m)
        assert pyo.value(m.sub.y) == pytest.approx(0.0)
