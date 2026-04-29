"""
Tests for the model validator (infeNMPC/check_model.py).

check_model runs three levels of validation on user-supplied model.py files:
  Level 1 — Static AST: valid Python, required functions present, no m.time
             declaration, both functions return m.
  Level 2 — Import: module loads without error.
  Level 3 — Dynamic Pyomo: variables_initialize and equations_write run against
             a real ConcreteModel; required components (MV_index, CV_index,
             setpoints, initial_values, DerivativeVar) are present.

These tests verify that:
  - The canonical CSTR model passes all three levels with no failures.
  - Intentionally broken models fail at the correct level with FAIL status.
  - Optional items (custom_objective, default_options) produce WARNs, not FAILs,
    so that models without economic objectives are not incorrectly rejected.
"""
import textwrap
from pathlib import Path

import pytest

from infeNMPC.check_model import validate

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_CSTR_MODEL = (
    Path(__file__).parent.parent
    / "examples"
    / "lyap_flag_examples"
    / "enmpc_cstr"
    / "small_cstr_model.py"
)


def _write_model(tmp_path, src: str) -> Path:
    """Write src to a temp model.py and return its Path."""
    p = tmp_path / "model.py"
    p.write_text(textwrap.dedent(src))
    return p


# ---------------------------------------------------------------------------
# Canonical CSTR model
# ---------------------------------------------------------------------------

class TestCSTRModelValidation:
    """The canonical CSTR model should pass all three validation levels."""

    def test_report_ok(self):
        """report.ok is True when there are no FAIL results."""
        report = validate(_CSTR_MODEL, dynamic=True)
        assert report.ok, [r for r in report.results if r.status == "FAIL"]

    def test_no_failures(self):
        """Zero individual check failures."""
        report = validate(_CSTR_MODEL, dynamic=True)
        assert report.failed == 0

    def test_has_passes(self):
        """At least the required static checks all pass."""
        report = validate(_CSTR_MODEL, dynamic=True)
        assert report.passed >= 1

    def test_static_level_passes(self):
        """Level 1 static check: syntax + required functions + no m.time."""
        report = validate(_CSTR_MODEL, dynamic=False)
        assert report.failed == 0

    def test_dynamic_level_passes(self):
        """Level 3 dynamic check: live Pyomo execution produces required components."""
        report = validate(_CSTR_MODEL, dynamic=True)
        assert report.failed == 0


# ---------------------------------------------------------------------------
# Models that should fail Level 1 (static AST)
# ---------------------------------------------------------------------------

class TestStaticFailures:
    """Broken Python or missing structure triggers FAIL at Level 1."""

    def test_syntax_error_fails(self, tmp_path):
        """A file with invalid Python syntax must fail the static check."""
        p = _write_model(tmp_path, "def this is not valid python {{{")
        report = validate(p, dynamic=False)
        assert not report.ok

    def test_missing_variables_initialize_fails(self, tmp_path):
        """variables_initialize() is required; absence must be a FAIL."""
        src = """
            import pyomo.environ as pyo
            def equations_write(m):
                return m
        """
        p = _write_model(tmp_path, src)
        report = validate(p, dynamic=False)
        assert not report.ok

    def test_missing_equations_write_fails(self, tmp_path):
        """equations_write() is required; absence must be a FAIL."""
        src = """
            import pyomo.environ as pyo
            def variables_initialize(m):
                return m
        """
        p = _write_model(tmp_path, src)
        report = validate(p, dynamic=False)
        assert not report.ok

    def test_declares_m_time_fails(self, tmp_path):
        """
        The framework injects m.time before calling variables_initialize.
        If the model also declares m.time, it would conflict.  The static check
        must catch this and produce a FAIL before any IPOPT is invoked.
        """
        src = """
            import pyomo.environ as pyo
            from pyomo.dae import ContinuousSet, DerivativeVar
            def variables_initialize(m):
                m.time = ContinuousSet(bounds=(0, 1))  # WRONG
                m.MV_index = pyo.Set(initialize=["u"])
                m.CV_index = pyo.Set(initialize=["x"])
                m.x = pyo.Var(m.time)
                m.dxdt = DerivativeVar(m.x, wrt=m.time)
                m.setpoints = pyo.Param(m.CV_index, initialize={"x": 0})
                m.initial_values = pyo.Param(m.CV_index, initialize={"x": 0})
                return m
            def equations_write(m):
                return m
        """
        p = _write_model(tmp_path, src)
        report = validate(p, dynamic=False)
        assert not report.ok

    def test_variables_initialize_missing_return_is_warned(self, tmp_path):
        """
        variables_initialize must return m.  Without the return, Pyomo state
        built inside the function is invisible to the framework.  The validator
        issues a WARN (not FAIL) so that models can still be discovered even
        when the static check fires — the dynamic check will catch actual failures.
        """
        src = """
            import pyomo.environ as pyo
            def variables_initialize(m):
                m.MV_index = pyo.Set(initialize=[])
                m.CV_index = pyo.Set(initialize=[])
                # no return
            def equations_write(m):
                return m
        """
        p = _write_model(tmp_path, src)
        report = validate(p, dynamic=False)
        # No FAIL — static level issues a WARN for missing return
        assert report.failed == 0
        assert report.warned >= 1

    def test_equations_write_missing_return_is_warned(self, tmp_path):
        """equations_write missing return also produces a WARN, not FAIL."""
        src = """
            import pyomo.environ as pyo
            def variables_initialize(m):
                m.MV_index = pyo.Set(initialize=[])
                m.CV_index = pyo.Set(initialize=[])
                return m
            def equations_write(m):
                pass  # no return
        """
        p = _write_model(tmp_path, src)
        report = validate(p, dynamic=False)
        assert report.failed == 0
        assert report.warned >= 1


# ---------------------------------------------------------------------------
# Optional items produce WARNs, not FAILs
# ---------------------------------------------------------------------------

class TestOptionalWarnings:
    """
    custom_objective and default_options are optional — their absence should
    produce a WARN (so users know about the feature) but not a FAIL (so that
    tracking-only models aren't broken).
    """

    def test_missing_custom_objective_is_warn_not_fail(self, tmp_path):
        src = """
            import pyomo.environ as pyo
            from pyomo.dae import DerivativeVar
            def variables_initialize(m):
                m.MV_index = pyo.Set(initialize=["u"])
                m.CV_index = pyo.Set(initialize=["x"])
                m.x = pyo.Var(m.time)
                m.dxdt = DerivativeVar(m.x, wrt=m.time)
                m.u = pyo.Var(m.time, bounds=(0, 1))
                m.setpoints = pyo.Param(m.CV_index, initialize={"x": 0.5})
                m.initial_values = pyo.Param(m.CV_index, initialize={"x": 0.0})
                return m
            def equations_write(m):
                m.ode = pyo.Constraint(m.time, rule=lambda m, t: m.dxdt[t] == -m.x[t])
                return m
        """
        p = _write_model(tmp_path, src)
        report = validate(p, dynamic=False)
        # No FAIL — just a WARN for the missing custom_objective
        assert report.failed == 0
