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
