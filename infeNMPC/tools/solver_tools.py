"""Solver utility functions: bound clipping and revive-run retry logic."""
import pyomo.environ as pyo
from pyomo.opt import TerminationCondition
from .debug_tools import _report_constraint_violations

_SKIP_REVIVE = frozenset({
    TerminationCondition.infeasible,
    TerminationCondition.maxIterations,
})


def _clip_to_bounds(model):
    """Project all variable values to [lb, ub].

    Handles tiny out-of-bound values introduced by bound_relax_factor != 0 in
    the warm solver. Prevents downstream Pyomo set_value() failures on
    NonNegativeReals slacks when loading data into subsequent models.
    """
    for vd in model.component_data_objects(pyo.Var, active=True, descend_into=True):
        val = pyo.value(vd)
        if val is None:
            continue
        if vd.lb is not None and val < vd.lb:
            vd.set_value(vd.lb)
        if vd.ub is not None and val > vd.ub:
            vd.set_value(vd.ub)


def _attempt_revive(controller, i, options, caught_exc):
    """Retry a failed controller solve up to options.revive_run times.

    If revive_run is 0 or the termination condition is not retryable (infeasible
    or maxIterations), re-raises caught_exc immediately.  Otherwise retries using
    controller.revive_solve() and prints status messages.  On success, returns
    normally (caller should then call _clip_to_bounds).  On exhaustion, re-raises
    the last RuntimeError from the final failed revive attempt.

    Parameters
    ----------
    controller : Controller
        The controller that just failed to solve.
    i : int
        Current MPC iteration index (for print messages).
    options : Options
        MPC options; reads revive_run and debug_flag.
    caught_exc : RuntimeError
        The exception raised by controller.solve().
    """
    tc = getattr(controller, 'last_tc', None)
    if options.revive_run <= 0 or tc in _SKIP_REVIVE:
        if options.debug_flag:
            _report_constraint_violations(
                controller._model, label=f"controller failure at iter {i}"
            )
        raise caught_exc

    for _attempt in range(options.revive_run):
        print(
            f"[revive_run] iter {i}: warm solve failed ({tc}), "
            f"retrying (attempt {_attempt + 1}/{options.revive_run})..."
        )
        try:
            controller.revive_solve()
            print(f"[revive_run] iter {i}: revival succeeded on attempt {_attempt + 1}.")
            return
        except RuntimeError:
            if _attempt == options.revive_run - 1:
                print(f"[revive_run] iter {i}: all {options.revive_run} attempt(s) exhausted.")
                if options.debug_flag:
                    _report_constraint_violations(
                        controller._model, label=f"revive failure at iter {i}"
                    )
                raise
            tc = getattr(controller, 'last_tc', tc)
