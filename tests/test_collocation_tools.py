"""
Tests for infeNMPC/tools/collocation_tools.py.

The collocation tools support the infinite-horizon time transformation.  The
key function is _compute_gamma_from_collocation, which auto-selects the
time-compression parameter γ so that the first LEGENDRE collocation point τ₁
satisfies τ₁ = tanh(γ · sampling_time).  This ensures the time-transformation
is calibrated to the actual collocation grid.

Also tested: _lagrange_coeffs_at_endpoint, which computes extrapolation
coefficients for evaluating algebraic CVs and MVs at τ=1 (a non-collocation
point in LEGENDRE discretization).
"""
import math
import pytest

from infeNMPC.tools.collocation_tools import (
    _compute_gamma_from_collocation,
    _get_last_fe_pts,
    _lagrange_coeffs_at_endpoint,
)


# ---------------------------------------------------------------------------
# _compute_gamma_from_collocation
# ---------------------------------------------------------------------------

class TestComputeGamma:
    """gamma is auto-computed from the collocation grid."""

    def test_returns_positive_float(self):
        """γ must be strictly positive for the time transform to map [0,∞)→[0,1)."""
        gamma = _compute_gamma_from_collocation(nfe=3, ncp=3, sampling_time=1.0)
        assert isinstance(gamma, float)
        assert gamma > 0

    def test_formula_atanh_tau1_over_sampling_time(self):
        """
        By design, γ = atanh(τ₁) / sampling_time.  Verify by rebuilding
        the same dummy model and checking tanh(γ · Ts) == τ₁.
        """
        import pyomo.environ as pyo
        from pyomo.dae import ContinuousSet, DerivativeVar

        nfe, ncp, Ts = 3, 3, 1.0
        gamma = _compute_gamma_from_collocation(nfe, ncp, Ts)

        dummy = pyo.ConcreteModel()
        dummy.time = ContinuousSet(bounds=(0, 1))
        dummy.x = pyo.Var(dummy.time)
        dummy.dxdt = DerivativeVar(dummy.x, wrt=dummy.time)
        pyo.TransformationFactory("dae.collocation").apply_to(
            dummy, ncp=ncp, nfe=nfe, wrt=dummy.time, scheme="LAGRANGE-LEGENDRE"
        )
        tau_1 = sorted(t for t in dummy.time if 0 < t < 1)[0]
        assert math.isclose(math.tanh(gamma * Ts), tau_1, rel_tol=1e-10)

    def test_scales_with_sampling_time(self):
        """
        Doubling the sampling time halves γ, because τ₁ is fixed by the
        collocation grid but γ = atanh(τ₁) / Ts.
        """
        g1 = _compute_gamma_from_collocation(3, 3, 1.0)
        g2 = _compute_gamma_from_collocation(3, 3, 2.0)
        assert math.isclose(g1 / g2, 2.0, rel_tol=1e-10)

    def test_different_nfe_gives_different_gamma(self):
        """
        More finite elements changes the collocation-point positions and
        therefore the auto-selected γ.
        """
        g1 = _compute_gamma_from_collocation(1, 3, 1.0)
        g2 = _compute_gamma_from_collocation(5, 3, 1.0)
        assert not math.isclose(g1, g2, rel_tol=1e-6)

    def test_nfe1_ncp1_produces_finite_gamma(self):
        """Even the coarsest grid (1 element, 1 point) gives a finite γ."""
        gamma = _compute_gamma_from_collocation(nfe=1, ncp=1, sampling_time=1.0)
        assert math.isfinite(gamma)
        assert gamma > 0


# ---------------------------------------------------------------------------
# _lagrange_coeffs_at_endpoint
# ---------------------------------------------------------------------------

class TestLagrangeCoeffs:
    """
    Lagrange extrapolation coefficients are used to evaluate algebraic CVs and
    MVs at τ=1, which is not a LEGENDRE collocation point.  The polynomial
    passes exactly through the collocation values, so:
      - coefficients sum to 1 (partition of unity at t_end=1 when the basis
        polynomial spans the same domain).
      - a single collocation point always maps to coefficient [1.0].
    """

    def test_single_point_coefficient_is_one(self):
        """With one collocation point, the only basis polynomial is L(t)=1."""
        coeffs = _lagrange_coeffs_at_endpoint([0.5], t_end=1.0)
        assert len(coeffs) == 1
        assert math.isclose(coeffs[0], 1.0, rel_tol=1e-12)

    def test_two_points_sum_to_one(self):
        """Lagrange basis polynomials form a partition of unity."""
        coeffs = _lagrange_coeffs_at_endpoint([0.2, 0.8], t_end=1.0)
        assert math.isclose(sum(coeffs), 1.0, rel_tol=1e-12)

    def test_three_points_sum_to_one(self):
        coeffs = _lagrange_coeffs_at_endpoint([0.1, 0.5, 0.9], t_end=1.0)
        assert math.isclose(sum(coeffs), 1.0, rel_tol=1e-12)

    def test_exact_recovery_at_collocation_point(self):
        """
        If t_end equals one of the collocation points, only that basis
        polynomial is non-zero (Kronecker-delta property).
        """
        pts = [0.2, 0.6]
        c = _lagrange_coeffs_at_endpoint(pts, t_end=0.2)
        assert math.isclose(c[0], 1.0, rel_tol=1e-12)
        assert math.isclose(c[1], 0.0, abs_tol=1e-12)

    def test_returns_correct_length(self):
        """Output length matches the number of input collocation points."""
        for n in [1, 2, 3, 5]:
            pts = [i / (n + 1) for i in range(1, n + 1)]
            coeffs = _lagrange_coeffs_at_endpoint(pts, t_end=1.0)
            assert len(coeffs) == n
