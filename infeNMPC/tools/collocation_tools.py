"""
Utilities for working with collocation discretizations in pyomo.dae.
"""
import math
import pyomo.environ as pyo
from pyomo.dae import ContinuousSet, DerivativeVar


def _compute_gamma_from_collocation(nfe, ncp, sampling_time):
    """
    Compute gamma so that the first infinite-block collocation point satisfies
    ``tau_1 = tanh(gamma * sampling_time)``.

    Builds a minimal dummy model, discretizes it with the same
    LAGRANGE-LEGENDRE scheme used for the infinite block, reads off the first
    non-boundary time point ``tau_1``, and returns
    ``atanh(tau_1) / sampling_time``.

    Parameters
    ----------
    nfe : int
        Number of finite elements for the infinite-horizon block.
    ncp : int
        Number of collocation points per finite element.
    sampling_time : float

    Returns
    -------
    float
    """
    dummy = pyo.ConcreteModel()
    dummy.time = ContinuousSet(bounds=(0, 1))
    dummy.x = pyo.Var(dummy.time)
    dummy.dxdt = DerivativeVar(dummy.x, wrt=dummy.time)
    discretizer = pyo.TransformationFactory('dae.collocation')
    discretizer.apply_to(
        dummy, ncp=ncp, nfe=nfe, wrt=dummy.time, scheme='LAGRANGE-LEGENDRE'
    )
    tau_1 = sorted(t for t in dummy.time if 0 < t < 1)[0]
    return math.atanh(tau_1) / sampling_time


def _get_last_fe_pts(block):
    """
    Return a sorted list of the collocation points that lie strictly inside
    the last finite element of *block*.

    Parameters
    ----------
    block : pyo.Block or pyo.ConcreteModel
        A discretized block with a ``time`` ContinuousSet.

    Returns
    -------
    list of float
    """
    fe_pts = list(block.time.get_finite_elements())
    t_start = fe_pts[-2]   # left boundary of the last FE
    t_end   = fe_pts[-1]   # right boundary (τ=1 for the infinite block)
    return sorted(t for t in block.time if t_start < t < t_end)


def _lagrange_coeffs_at_endpoint(colloc_pts, t_end):
    """
    Compute Lagrange basis polynomial coefficients evaluated at *t_end*.

    Given *n* collocation points ``τ_0, …, τ_{n-1}``, the *j*-th coefficient
    is the Lagrange basis polynomial ``L_j`` evaluated at *t_end*:

    .. math::

        L_j(t_{\\text{end}}) = \\prod_{k \\neq j}
            \\frac{t_{\\text{end}} - \\tau_k}{\\tau_j - \\tau_k}

    These are pure numerical scalars (independent of variable values), so they
    can be precomputed once and reused when building Pyomo expressions.

    Parameters
    ----------
    colloc_pts : list of float
        Collocation points in the last finite element (strictly interior to
        the element, i.e. excluding the element boundaries).
    t_end : float
        Target evaluation point (typically 1.0 for the infinite-horizon block).

    Returns
    -------
    list of float
        Length equals ``len(colloc_pts)``.
    """
    coeffs = []
    for j, tj in enumerate(colloc_pts):
        c = 1.0
        for k, tk in enumerate(colloc_pts):
            if k != j:
                c *= (t_end - tk) / (tj - tk)
        coeffs.append(c)
    return coeffs
