"""
Plant model for closed-loop NMPC simulations.
"""
import pyomo.environ as pyo
from .make_model import _make_finite_horizon_model, _ipopt_solver, _check_optimal
from .infNMPC_options import Options


class Plant:
    """
    A single-step dynamic simulation model used as the closed-loop plant.

    Wraps a Pyomo ``ConcreteModel`` built with finite-horizon settings
    (``nfe_finite=1``). Manipulated variables are fixed so the plant just
    integrates forward one sampling interval given the current MV values.

    Attribute access on a ``Plant`` instance falls through to the underlying
    Pyomo model, so all existing ``plant.CV_index``, ``plant.interface``,
    ``plant.time``, etc. accesses work unchanged.

    Parameters
    ----------
    options : Options
        Simulation configuration.  ``nfe_finite`` and ``infinite_horizon`` are
        overridden internally; all other fields are respected.

    Attributes
    ----------
    options : Options
        The options used to build this plant.
    """

    def __init__(self, options: Options):
        self.options = options
        plant_options = options.copy(nfe_finite=1, infinite_horizon=False)

        m = pyo.ConcreteModel()
        m = _make_finite_horizon_model(m, plant_options)

        print('Generating Plant')

        for var_name in m.MV_index:
            var = getattr(m, var_name)
            var.fix()
            if options.ncp_finite > 1:
                getattr(m, f"{var_name}_interpolation_constraints").deactivate()

        m.obj = pyo.Objective(expr=1)

        self._model = m
        self._solver = _ipopt_solver()

        print('Plant Initial Solve')
        self.solve()

    def solve(self):
        """Solve the plant model in place."""
        results = self._solver.solve(self._model, tee=self.options.tee_flag)
        _check_optimal(results, "plant")

    def __getattr__(self, name: str):
        # Delegate unknown attribute lookups to the underlying Pyomo model.
        # __getattr__ is only called when normal lookup fails, so 'options',
        # '_model', '_solver', and 'solve' are all found first.
        return getattr(self._model, name)
