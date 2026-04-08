"""
Debugging utilities for infeNMPC Pyomo models.
"""
import pyomo.environ as pyo


def _report_constraint_violations(model, n=10, tol=1e-6, label=""):
    """
    Print the *n* most violated active constraints in a Pyomo model.

    Traverses all active ``Constraint`` components recursively (including
    nested Blocks).  For each constraint data object the violation is:

    * equality  (lower == upper): ``|body - lower|``
    * inequality (lower <= body <= upper):
      ``max(0, lower - body, body - upper)``

    Only violations larger than *tol* are reported.

    Parameters
    ----------
    model : pyo.ConcreteModel or pyo.Block
        The model to inspect.
    n : int
        Maximum number of violations to print.
    tol : float
        Minimum violation magnitude to include in the report.
    label : str
        Optional label appended to the header line.
    """
    violations = []

    for con in model.component_objects(pyo.Constraint, active=True, descend_into=True):
        for idx in con:
            condata = con[idx]
            try:
                body_val = pyo.value(condata.body, exception=False)
            except Exception:
                continue
            if body_val is None:
                continue

            lower = pyo.value(condata.lower) if condata.lower is not None else None
            upper = pyo.value(condata.upper) if condata.upper is not None else None

            if lower is not None and upper is not None and lower == upper:
                # equality constraint
                viol = abs(body_val - lower)
            else:
                viol = 0.0
                if lower is not None:
                    viol = max(viol, lower - body_val)
                if upper is not None:
                    viol = max(viol, body_val - upper)

            if viol > tol:
                violations.append((viol, con.name, idx))

    violations.sort(key=lambda x: -x[0])

    header = "Constraint violation report"
    if label:
        header += f"  [{label}]"
    print(f"\n{'='*60}")
    print(header)
    print(f"{'='*60}")
    if not violations:
        print(f"  No violations above tol={tol:.1e}")
    else:
        print(f"  {len(violations)} violation(s) found; showing top {min(n, len(violations))}:")
        print(f"  {'Violation':>14}  Constraint")
        print(f"  {'-'*14}  {'-'*40}")
        for viol, name, idx in violations[:n]:
            idx_str = str(idx) if not isinstance(idx, tuple) else ", ".join(str(i) for i in idx)
            print(f"  {viol:>14.6g}  {name}[{idx_str}]")
    print(f"{'='*60}\n")
