"""
check_model.py — Validator for infeNMPC model.py files.

Runs three levels of checks:
  1. Static  — AST parse, required function defs, return statements, m.time guard
  2. Import  — module loads cleanly
  3. Dynamic — executes variables_initialize and equations_write against a real
               Pyomo ConcreteModel and inspects the resulting components

Usage (CLI entry point):
    check-model path/to/model.py
    check-model path/to/model.py --no-dynamic
"""
from __future__ import annotations

import ast
import importlib.util
import sys
import textwrap
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


# ──────────────────────────────────────────────────────────────────────────────
# Result types
# ──────────────────────────────────────────────────────────────────────────────

PASS = "PASS"
FAIL = "FAIL"
WARN = "WARN"
SKIP = "SKIP"


@dataclass
class CheckResult:
    status: str          # PASS | FAIL | WARN | SKIP
    name: str
    detail: str = ""


@dataclass
class ValidationReport:
    results: list[CheckResult] = field(default_factory=list)

    def add(self, status: str, name: str, detail: str = ""):
        self.results.append(CheckResult(status, name, detail))

    @property
    def passed(self) -> int:
        return sum(1 for r in self.results if r.status == PASS)

    @property
    def failed(self) -> int:
        return sum(1 for r in self.results if r.status == FAIL)

    @property
    def warned(self) -> int:
        return sum(1 for r in self.results if r.status == WARN)

    @property
    def skipped(self) -> int:
        return sum(1 for r in self.results if r.status == SKIP)

    @property
    def ok(self) -> bool:
        return self.failed == 0


# ──────────────────────────────────────────────────────────────────────────────
# ANSI colours (disabled on non-TTY)
# ──────────────────────────────────────────────────────────────────────────────

def _colour(code: str, text: str) -> str:
    if sys.stdout.isatty():
        return f"\033[{code}m{text}\033[0m"
    return text

def _green(t):  return _colour("32", t)
def _red(t):    return _colour("31", t)
def _yellow(t): return _colour("33", t)
def _dim(t):    return _colour("2",  t)


def _badge(status: str) -> str:
    if status == PASS: return _green("[PASS]")
    if status == FAIL: return _red  ("[FAIL]")
    if status == WARN: return _yellow("[WARN]")
    return _dim("[SKIP]")


# ──────────────────────────────────────────────────────────────────────────────
# Level 1: Static AST checks
# ──────────────────────────────────────────────────────────────────────────────

def _ast_checks(source: str, report: ValidationReport) -> Optional[ast.Module]:
    """Parse and run AST-level checks. Returns the tree, or None on parse error."""
    # 1a. Valid Python syntax
    try:
        tree = ast.parse(source)
        report.add(PASS, "Valid Python syntax")
    except SyntaxError as exc:
        report.add(FAIL, "Valid Python syntax", f"SyntaxError: {exc}")
        return None

    # Collect top-level function definitions
    top_level_fns: dict[str, ast.FunctionDef] = {
        node.name: node
        for node in ast.walk(tree)
        if isinstance(node, ast.FunctionDef) and _is_top_level(node, tree)
    }

    # 1b. Required functions present
    for fn_name in ("variables_initialize", "equations_write"):
        if fn_name in top_level_fns:
            report.add(PASS, f"{fn_name}() defined")
        else:
            report.add(FAIL, f"{fn_name}() defined", "Function not found in module")

    # 1c. Optional functions (warn if absent)
    for fn_name in ("custom_objective", "default_options"):
        if fn_name in top_level_fns:
            report.add(PASS, f"{fn_name}() defined (optional)")
        else:
            report.add(WARN, f"{fn_name}() defined (optional)", "Not present — only needed if custom_objective=True or model-level defaults required")

    # 1d. variables_initialize must NOT declare m.time
    if "variables_initialize" in top_level_fns:
        fn_node = top_level_fns["variables_initialize"]
        if _declares_m_time(fn_node):
            report.add(
                FAIL,
                "variables_initialize() does not declare m.time",
                "m.time is injected by the framework — declaring it here will cause a conflict",
            )
        else:
            report.add(PASS, "variables_initialize() does not declare m.time")

    # 1e. Both required functions return m
    for fn_name in ("variables_initialize", "equations_write"):
        if fn_name not in top_level_fns:
            report.add(SKIP, f"{fn_name}() returns m", "Function not found")
            continue
        fn_node = top_level_fns[fn_name]
        if _has_return_m(fn_node):
            report.add(PASS, f"{fn_name}() returns m")
        else:
            report.add(WARN, f"{fn_name}() returns m", "No 'return m' found — framework expects both functions to return the model")

    return tree


def _is_top_level(node: ast.FunctionDef, tree: ast.Module) -> bool:
    """True if node is a direct child of the module body (not nested)."""
    return any(child is node for child in ast.iter_child_nodes(tree))


def _declares_m_time(fn_node: ast.FunctionDef) -> bool:
    """
    Return True if the function body contains any assignment that looks like
    m.time = ... or m.time: ... = ...
    """
    for node in ast.walk(fn_node):
        # m.time = <value>  →  Assign with target being Attribute(attr='time')
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if _is_attr(target, "time"):
                    return True
        # m.time: type = value  →  AnnAssign
        if isinstance(node, ast.AnnAssign):
            if _is_attr(node.target, "time"):
                return True
    return False


def _is_attr(node: ast.expr, attr: str) -> bool:
    return isinstance(node, ast.Attribute) and node.attr == attr


def _has_return_m(fn_node: ast.FunctionDef) -> bool:
    """True if there is at least one 'return m' statement in the function."""
    for node in ast.walk(fn_node):
        if isinstance(node, ast.Return) and node.value is not None:
            if isinstance(node.value, ast.Name) and node.value.id == "m":
                return True
    return False


# ──────────────────────────────────────────────────────────────────────────────
# Level 2: Import check
# ──────────────────────────────────────────────────────────────────────────────

def _import_check(path: Path, report: ValidationReport):
    """Try to import the module. Returns the live module or None."""
    sys.path.insert(0, str(path.parent))
    spec = importlib.util.spec_from_file_location(path.stem, path)
    mod = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(mod)
        report.add(PASS, "Module imports without error")
        return mod
    except Exception as exc:
        report.add(FAIL, "Module imports without error", _truncate(repr(exc)))
        return None


# ──────────────────────────────────────────────────────────────────────────────
# Level 3: Dynamic Pyomo checks
# ──────────────────────────────────────────────────────────────────────────────

_PYOMO_IMPORT_MSG = "Pyomo not installed or import failed — skipping dynamic checks"


def _dynamic_checks(mod, report: ValidationReport):
    """Run variables_initialize and equations_write on a real Pyomo model."""
    # Try importing Pyomo
    try:
        import pyomo.environ as pyo
        from pyomo.dae import ContinuousSet, DerivativeVar
        from pyomo.environ import TransformationFactory
    except ImportError:
        report.add(SKIP, "Dynamic checks", _PYOMO_IMPORT_MSG)
        return

    # ── variables_initialize ──────────────────────────────────────────────────
    vi_fn = getattr(mod, "variables_initialize", None)
    if vi_fn is None:
        report.add(SKIP, "variables_initialize() runs on ConcreteModel", "Function not found")
        m_after_vi = None
    else:
        m = pyo.ConcreteModel()
        m.time = ContinuousSet(bounds=(0, 2))
        try:
            result = vi_fn(m)
            report.add(PASS, "variables_initialize() runs on ConcreteModel")
            m_after_vi = result if result is not None else m
        except Exception as exc:
            report.add(FAIL, "variables_initialize() runs on ConcreteModel", _truncate(repr(exc)))
            m_after_vi = None

    if m_after_vi is not None:
        m = m_after_vi

        # Required component checks
        for attr, label in [
            ("MV_index", "m.MV_index declared"),
            ("CV_index", "m.CV_index declared"),
            ("setpoints", "m.setpoints declared"),
            ("initial_values", "m.initial_values declared"),
        ]:
            if hasattr(m, attr):
                report.add(PASS, label)
            else:
                report.add(FAIL, label, f"m.{attr} not found after variables_initialize()")

        # At least one DerivativeVar — use component_objects(DerivativeVar) directly
        try:
            deriv_vars = list(m.component_objects(DerivativeVar))
        except Exception:
            # Fallback: walk all components and check type name
            deriv_vars = [
                comp for comp in m.component_objects(pyo.Var)
                if type(comp).__name__ == "DerivativeVar"
            ]
        if deriv_vars:
            names = ", ".join(c.name for c in deriv_vars)
            report.add(PASS, "At least one DerivativeVar declared", names)
        else:
            report.add(WARN, "At least one DerivativeVar declared",
                       "No DerivativeVar found — algebraic-only models are uncommon but allowed")

        # MV_index and CV_index non-empty
        for attr in ("MV_index", "CV_index"):
            if hasattr(m, attr):
                idx = getattr(m, attr)
                try:
                    members = list(idx)
                    if members:
                        report.add(PASS, f"m.{attr} non-empty", f"Members: {members}")
                    else:
                        report.add(WARN, f"m.{attr} non-empty", f"m.{attr} is empty")
                except Exception:
                    pass

        # ── equations_write ───────────────────────────────────────────────────
        ew_fn = getattr(mod, "equations_write", None)
        if ew_fn is None:
            report.add(SKIP, "equations_write() runs after discretization", "Function not found")
        else:
            # Discretize with simple LAGRANGE-RADAU (same as finite block)
            try:
                discretizer = TransformationFactory("dae.collocation")
                discretizer.apply_to(m, wrt=m.time, nfe=2, ncp=1, scheme="LAGRANGE-RADAU")
            except Exception as exc:
                report.add(
                    WARN,
                    "equations_write() runs after discretization",
                    f"Could not discretize model for dynamic check: {_truncate(repr(exc))}",
                )
                return

            try:
                ew_fn(m)
                report.add(PASS, "equations_write() runs after discretization")
            except Exception as exc:
                report.add(FAIL, "equations_write() runs after discretization", _truncate(repr(exc)))


# ──────────────────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────────────────

def _truncate(s: str, limit: int = 200) -> str:
    return s if len(s) <= limit else s[:limit] + "…"


# ──────────────────────────────────────────────────────────────────────────────
# Public API
# ──────────────────────────────────────────────────────────────────────────────

def validate(path: Path, dynamic: bool = True) -> ValidationReport:
    """
    Validate a model.py file and return a ValidationReport.

    Parameters
    ----------
    path:
        Absolute or relative path to the model.py file.
    dynamic:
        If True, run dynamic Pyomo checks (requires Pyomo to be installed).
    """
    report = ValidationReport()
    path = Path(path).resolve()

    # File existence
    if not path.exists():
        report.add(FAIL, "File exists", str(path))
        return report
    report.add(PASS, "File exists")

    source = path.read_text(encoding="utf-8")

    # Level 1
    _ast_checks(source, report)

    # Level 2
    mod = _import_check(path, report)

    # Level 3
    if mod is not None and dynamic:
        _dynamic_checks(mod, report)
    elif not dynamic:
        report.add(SKIP, "Dynamic checks", "Skipped via --no-dynamic")

    return report


# ──────────────────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────────────────

def _print_report(report: ValidationReport, path: Path):
    print(f"\nValidating {path}\n" + "─" * 60)
    for r in report.results:
        line = f"  {_badge(r.status)}  {r.name}"
        print(line)
        if r.detail:
            wrapped = textwrap.fill(r.detail, width=72, initial_indent="          ", subsequent_indent="          ")
            print(_dim(wrapped))

    total = len(report.results)
    summary_parts = [f"{report.passed}/{total} checks passed"]
    if report.warned:
        summary_parts.append(f"{report.warned} warning{'s' if report.warned > 1 else ''}")
    if report.skipped:
        summary_parts.append(f"{report.skipped} skipped")

    print("\n" + "─" * 60)
    if report.ok:
        print(_green("✓ " + ", ".join(summary_parts)))
    else:
        print(_red(f"✗ {report.failed} check{'s' if report.failed > 1 else ''} failed  ·  " + ", ".join(summary_parts)))
    print()


def main():
    import argparse

    parser = argparse.ArgumentParser(
        prog="check-model",
        description="Validate an infeNMPC model.py file.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""\
            Checks performed
            ────────────────
            Static (AST):
              • Valid Python syntax
              • variables_initialize() and equations_write() are defined
              • custom_objective() and default_options() detected (optional)
              • variables_initialize() does not declare m.time (injected by framework)
              • Both required functions return m

            Import:
              • Module loads without error

            Dynamic (Pyomo):
              • variables_initialize() runs on a ConcreteModel with m.time set
              • m.MV_index, m.CV_index, m.setpoints, m.initial_values are declared
              • At least one DerivativeVar is present
              • m.MV_index and m.CV_index are non-empty
              • equations_write() runs after collocation discretization
        """),
    )
    parser.add_argument(
        "model",
        metavar="MODEL_FILE",
        help="Path to a model.py file (e.g. examples/enmpc_cstr/model.py)",
    )
    parser.add_argument(
        "--no-dynamic",
        action="store_true",
        help="Skip the dynamic Pyomo checks (faster; no Pyomo required)",
    )
    args = parser.parse_args()

    path = Path(args.model)
    report = validate(path, dynamic=not args.no_dynamic)
    _print_report(report, path)
    sys.exit(0 if report.ok else 1)


if __name__ == "__main__":
    main()
