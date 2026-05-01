"""
Debugging utilities for infeNMPC Pyomo models.
"""
import pyomo.environ as pyo

# Threshold separating "truly optimal" IPOPT residuals (~1e-10) from
# "acceptable" residuals (~1e-7).  Any equality-constraint violation above
# this value at a nominally-successful solve indicates IPOPT exited via the
# acceptable-level criterion rather than true KKT optimality.
_IPOPT_OPTIMAL_TOL = 5e-8


def _report_acceptable_termination(model, results=None, label=""):
    """
    Detect and report constraints that prevented true optimal convergence.

    When IPOPT exits via the acceptable-level criterion (inf_pr ~1e-7)
    instead of true KKT optimality (inf_pr ~1e-10), this function prints
    a clear warning and lists the equality-constraint violations responsible.

    Detection is attempted in two ways, in order:
    1. The IPOPT solver-exit message is inspected for the word "acceptable".
    2. Fallback: any active constraint with residual > ``_IPOPT_OPTIMAL_TOL``
       (5e-8) is taken as evidence of acceptable termination.

    Parameters
    ----------
    model : pyo.ConcreteModel or pyo.Block
        The model that was just solved.
    results : solver results object, optional
        The object returned by ``solver.solve()``.  When provided, the IPOPT
        exit message is checked first.  When absent or when the message
        attribute is unavailable, falls back to the constraint-violation check.
    label : str
        Optional label appended to the report header.

    Returns
    -------
    bool
        ``True`` if acceptable termination was detected and reported,
        ``False`` if the solve appears to have reached true optimality.
    """
    is_acceptable = False

    # ---- Primary detection: IPOPT solver exit message ----
    if results is not None:
        solver_obj = getattr(results, 'solver', None)
        if solver_obj is not None:
            for attr in ('message', 'termination_message'):
                try:
                    msg = getattr(solver_obj, attr, None) or ''
                    if 'acceptable' in str(msg).lower():
                        is_acceptable = True
                        break
                except Exception:
                    pass

    # ---- Fallback detection: constraint violations above optimal tolerance ----
    if not is_acceptable:
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
                    viol = abs(body_val - lower)
                else:
                    viol = 0.0
                    if lower is not None:
                        viol = max(viol, lower - body_val)
                    if upper is not None:
                        viol = max(viol, body_val - upper)
                if viol > _IPOPT_OPTIMAL_TOL:
                    is_acceptable = True
                    break
            if is_acceptable:
                break

    if is_acceptable:
        suffix = f"  [{label}]" if label else ""
        print(f"\n{'!'*60}")
        print(f"  Acceptable termination{suffix}")
        print(f"  Constraints preventing optimal convergence (residual > {_IPOPT_OPTIMAL_TOL:.0e}):")
        print(f"{'!'*60}")
        _report_constraint_violations(model, n=10, tol=_IPOPT_OPTIMAL_TOL, label=label)

    return is_acceptable


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


def _report_solver_diagnostics(controller, options, iteration, eps_tol=1e-4):
    """
    Print per-iteration diagnostics: spec slack values, terminal constraint
    status, and Lyapunov constraint status.  Called when ``debug_flag=True``.

    Parameters
    ----------
    controller : Controller
        The controller after a successful solve.
    options : Options
        Simulation configuration.
    iteration : int
        Current MPC iteration index (0-based).
    eps_tol : float
        Threshold for flagging spec slacks as violated (default 1e-4).
    """
    block = controller.finite_block if options.infinite_horizon else controller._model
    t0 = block.time.first()

    print(f"\n{'─'*60}")
    print(f"  MPC diagnostics  [iter {iteration}]")
    print(f"{'─'*60}")

    def _max_of_var(sv):
        """Max value across all indices, excluding t=t0 in the time dimension."""
        vals = []
        for idx in sv.index_set():
            time_val = idx[-1] if isinstance(idx, tuple) else idx
            if time_val == t0:
                continue
            v = pyo.value(sv[idx], exception=False)
            if v is not None:
                vals.append(v)
        return max(vals) if vals else 0.0

    # ---- Spec slacks (product-purity specification violations) ----
    if hasattr(block, 'spec_slack_index'):
        max_vals = {}
        for sv_name in block.spec_slack_index:
            sv = getattr(block, sv_name, None)
            if sv is not None and isinstance(sv, pyo.Var):
                max_vals[sv_name] = _max_of_var(sv)

        if not max_vals:
            print(f"  spec slacks: (none found)")
        elif all(v <= eps_tol for v in max_vals.values()):
            detail = "  ".join(f"{n}={v:.2e}" for n, v in max_vals.items())
            print(f"  spec slacks OK (all <= {eps_tol:.0e}):  {detail}")
        else:
            for sv_name, max_val in max_vals.items():
                flag = "  ** VIOLATED **" if max_val > eps_tol else "  ok"
                print(f"  spec {sv_name}: max={max_val:.3e}{flag}")
    else:
        print(f"  spec slacks: no spec_slack_index on block")

    # ---- Physical slacks (slack_index) ----
    if hasattr(block, 'slack_index'):
        max_vals = {}
        for sv_name in block.slack_index:
            sv = getattr(block, sv_name, None)
            if sv is not None and isinstance(sv, pyo.Var):
                max_vals[sv_name] = _max_of_var(sv)

        if not max_vals:
            print(f"  phys slacks: (none found)")
        elif all(v <= eps_tol for v in max_vals.values()):
            detail = "  ".join(f"{n}={v:.2e}" for n, v in max_vals.items())
            print(f"  phys  slacks OK (all <= {eps_tol:.0e}):  {detail}")
        else:
            for sv_name, max_val in max_vals.items():
                flag = "  ** VIOLATED **" if max_val > eps_tol else "  ok"
                print(f"  phys  {sv_name}: max={max_val:.3e}{flag}")

    # ---- Terminal constraint ----
    tct = options.terminal_constraint_type
    if tct == 'soft':
        if options.infinite_horizon:
            ib = getattr(controller, 'infinite_block', None)
            penalty = (pyo.value(ib.infinite_terminal_soft_penalty)
                       if ib is not None and hasattr(ib, 'infinite_terminal_soft_penalty')
                       else None)
        else:
            m = controller._model
            penalty = (pyo.value(m.finite_terminal_soft_penalty)
                       if hasattr(m, 'finite_terminal_soft_penalty')
                       else None)
        if penalty is not None:
            flag = "OK" if penalty <= eps_tol else "** VIOLATED **"
            print(f"  terminal (soft): unweighted penalty = {penalty:.3e}  [{flag}]")
        else:
            print(f"  terminal (soft): penalty expression not found")
    elif tct == 'hard':
        print(f"  terminal (hard): enforced as equality at optimal")
    else:
        print(f"  terminal: not active (type='{tct}')")

    # ---- Lyapunov constraint ----
    if options.lyap_flag:
        lyap_ct = getattr(options, 'lyap_constraint_type', 'hard')
        if lyap_ct == 'soft':
            m = controller._model
            if hasattr(m, 'lyap_slack'):
                lyap_val = pyo.value(m.lyap_slack)
                flag = "OK" if lyap_val <= eps_tol else "** VIOLATED **"
                print(f"  lyapunov (soft): slack = {lyap_val:.3e}  [{flag}]")
            else:
                print(f"  lyapunov (soft): lyap_slack not found on model")
        elif lyap_ct == 'hard':
            print(f"  lyapunov (hard): decrease constraint enforced at optimal")
        else:
            print(f"  lyapunov: infrastructure built, constraint not active (type='none')")
    else:
        print(f"  lyapunov: not active (lyap_flag=False)")

    print(f"{'─'*60}\n")
