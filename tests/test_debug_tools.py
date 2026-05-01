"""
Tests for infeNMPC/tools/debug_tools.py — specifically _report_acceptable_termination.

_report_acceptable_termination detects when an IPOPT solve terminated at
acceptable quality (inf_pr ~1e-7) rather than true KKT optimality
(inf_pr ~1e-10) and prints the violated constraints responsible.

Detection uses two mechanisms in priority order:
  1. IPOPT solver exit message contains the word "acceptable".
  2. Fallback: any active constraint has residual > _IPOPT_OPTIMAL_TOL (5e-8).

These tests use hand-built Pyomo models with known constraint values so that
IPOPT is not required.  The model variables are set to specific values to
simulate what the model state looks like after an acceptable-level IPOPT exit.
"""
import pytest
import pyomo.environ as pyo

from infeNMPC.tools.debug_tools import (
    _report_acceptable_termination,
    _IPOPT_OPTIMAL_TOL,
)


# ---------------------------------------------------------------------------
# Helper: minimal mock of a Pyomo solver results object
# ---------------------------------------------------------------------------

class _MockSolverResult:
    """Minimal stand-in for the object returned by solver.solve()."""

    def __init__(self, message=""):
        self.solver = _MockSolverAttr(message=message)


class _MockSolverAttr:
    def __init__(self, message=""):
        self.message = message


# ---------------------------------------------------------------------------
# Primary detection: solver message
# ---------------------------------------------------------------------------

class TestAcceptableTerminationViaMessage:
    """_report_acceptable_termination detects via IPOPT exit message."""

    def _trivial_model(self):
        """A model with no constraints — detection must rely on message alone."""
        m = pyo.ConcreteModel()
        m.x = pyo.Var(initialize=0.0)
        return m

    def test_message_acceptable_triggers_report(self, capsys):
        """'Solved To Acceptable Level.' in message → detected."""
        m = self._trivial_model()
        results = _MockSolverResult(message="Solved To Acceptable Level.")
        detected = _report_acceptable_termination(m, results=results, label="unit-test")
        assert detected is True
        out = capsys.readouterr().out
        assert "Acceptable termination" in out

    def test_message_optimal_no_report(self, capsys):
        """'Optimal Solution Found.' in message → NOT detected."""
        m = self._trivial_model()
        results = _MockSolverResult(message="Optimal Solution Found.")
        detected = _report_acceptable_termination(m, results=results, label="unit-test")
        assert detected is False
        out = capsys.readouterr().out
        assert "Acceptable termination" not in out

    def test_message_case_insensitive(self, capsys):
        """Detection is case-insensitive."""
        m = self._trivial_model()
        results = _MockSolverResult(message="SOLVED TO ACCEPTABLE LEVEL.")
        detected = _report_acceptable_termination(m, results=results)
        assert detected is True

    def test_none_results_falls_through_to_violation_check(self, capsys):
        """When results=None, fallback violation check is used."""
        m = pyo.ConcreteModel()
        m.x = pyo.Var(initialize=1e-7)
        m.eq = pyo.Constraint(expr=m.x == 0.0)  # violation = 1e-7 > 5e-8
        detected = _report_acceptable_termination(m, results=None)
        assert detected is True


# ---------------------------------------------------------------------------
# Fallback detection: constraint violations
# ---------------------------------------------------------------------------

class TestAcceptableTerminationViaViolation:
    """_report_acceptable_termination fallback: constraint violation > _IPOPT_OPTIMAL_TOL."""

    def test_violation_above_threshold_detected(self, capsys):
        """Equality constraint with residual 1e-7 > 5e-8 is flagged."""
        m = pyo.ConcreteModel()
        m.x = pyo.Var(initialize=1e-7)
        m.eq = pyo.Constraint(expr=m.x == 0.0)
        detected = _report_acceptable_termination(m)
        assert detected is True
        out = capsys.readouterr().out
        assert "Acceptable termination" in out

    def test_violation_below_threshold_not_detected(self, capsys):
        """Residual 1e-10 (well below 5e-8) is treated as truly optimal."""
        m = pyo.ConcreteModel()
        m.x = pyo.Var(initialize=1e-10)
        m.eq = pyo.Constraint(expr=m.x == 0.0)
        detected = _report_acceptable_termination(m)
        assert detected is False
        out = capsys.readouterr().out
        assert "Acceptable termination" not in out

    def test_no_constraints_no_report(self, capsys):
        """A model with no constraints and results=None → not acceptable."""
        m = pyo.ConcreteModel()
        m.x = pyo.Var(initialize=1e-3)  # large value, but no constraints
        detected = _report_acceptable_termination(m)
        assert detected is False

    def test_violation_exactly_at_threshold_not_detected(self):
        """Residual exactly equal to _IPOPT_OPTIMAL_TOL is NOT flagged (strictly >)."""
        m = pyo.ConcreteModel()
        m.x = pyo.Var(initialize=_IPOPT_OPTIMAL_TOL)
        m.eq = pyo.Constraint(expr=m.x == 0.0)
        detected = _report_acceptable_termination(m)
        assert detected is False

    def test_inequality_constraint_violation_detected(self, capsys):
        """Upper-bound inequality violation above threshold is also flagged."""
        m = pyo.ConcreteModel()
        m.x = pyo.Var(initialize=1.0 + 1e-7)  # 1e-7 above upper bound of 1.0
        m.ineq = pyo.Constraint(expr=m.x <= 1.0)
        detected = _report_acceptable_termination(m)
        assert detected is True

    def test_nested_block_constraint_detected(self, capsys):
        """Violation inside a nested Block is found (descend_into=True)."""
        m = pyo.ConcreteModel()
        m.sub = pyo.Block()
        m.sub.y = pyo.Var(initialize=1e-7)
        m.sub.eq = pyo.Constraint(expr=m.sub.y == 0.0)
        detected = _report_acceptable_termination(m)
        assert detected is True

    def test_user_suggested_lhs_rhs_case(self, capsys):
        """
        Specific case the user described: LHS=1e-8, RHS=1e-7.
        Constraint: x == 1e-7, but x is set to 1e-8.
        Violation = |1e-8 - 1e-7| = 9e-8 > 5e-8 → acceptable detected.
        """
        m = pyo.ConcreteModel()
        m.x = pyo.Var(initialize=1e-8)      # LHS value
        m.eq = pyo.Constraint(expr=m.x == 1e-7)  # RHS = 1e-7
        detected = _report_acceptable_termination(m)
        assert detected is True
        out = capsys.readouterr().out
        assert "Acceptable termination" in out
        assert "eq" in out  # the violated constraint name appears in the report


# ---------------------------------------------------------------------------
# Output format
# ---------------------------------------------------------------------------

class TestOutputFormat:
    """The printed output matches the expected format."""

    def test_label_appears_in_output(self, capsys):
        """The label parameter appears in the output header."""
        m = pyo.ConcreteModel()
        m.x = pyo.Var(initialize=1e-7)
        m.eq = pyo.Constraint(expr=m.x == 0.0)
        _report_acceptable_termination(m, label="my-label")
        out = capsys.readouterr().out
        assert "my-label" in out

    def test_return_false_for_truly_optimal(self):
        """Returns False (not just falsy) for a truly optimal model."""
        m = pyo.ConcreteModel()
        m.x = pyo.Var(initialize=1e-12)
        m.eq = pyo.Constraint(expr=m.x == 0.0)
        assert _report_acceptable_termination(m) is False

    def test_return_true_for_acceptable(self):
        """Returns True (not just truthy) for an acceptable-level model."""
        m = pyo.ConcreteModel()
        m.x = pyo.Var(initialize=1e-7)
        m.eq = pyo.Constraint(expr=m.x == 0.0)
        assert _report_acceptable_termination(m) is True

    def test_threshold_value_is_correct(self):
        """_IPOPT_OPTIMAL_TOL is 5e-8 as documented."""
        assert _IPOPT_OPTIMAL_TOL == pytest.approx(5e-8)
