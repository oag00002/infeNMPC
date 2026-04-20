"""
Plant model for closed-loop NMPC simulations.
"""
import pyomo.environ as pyo
from .make_model import _make_finite_horizon_model, _ipopt_solver, _check_optimal
from .infNMPC_options import Options
from .tools.debug_tools import _report_constraint_violations


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
        Simulation configuration.  ``nfe_finite``, ``infinite_horizon``, and
        ``terminal_constraint_type`` are overridden internally (plant is always
        a plain forward simulation with no terminal constraints); all other
        fields are respected.

    Attributes
    ----------
    options : Options
        The options used to build this plant.
    """

    def __init__(self, options: Options):
        self.options = options
        plant_options = options.copy(nfe_finite=1, infinite_horizon=False,
                                     terminal_constraint_type='none',
                                     lyap_flag=False)

        m = pyo.ConcreteModel()
        m = _make_finite_horizon_model(m, plant_options)

        print('Generating Plant')

        for var_name in m.MV_index:
            var = getattr(m, var_name)
            var.fix()
            if options.ncp_finite > 1:
                getattr(m, f"{var_name}_interpolation_constraints").deactivate()

        # Deactivate interpolation constraints for slack variables (same as MVs).
        # Do NOT fix them here — DynamicModelInterface.load_data skips fixed
        # variables, so slacks must remain free when the warm-start is loaded in
        # run_MPC.py.  After load_data the caller must call sv_var.fix() so that
        # the plant constraints exactly match the controller's first FE.
        if hasattr(m, 'slack_index'):
            for sv_name in m.slack_index:
                sv_var = getattr(m, sv_name, None)
                if sv_var is not None and isinstance(sv_var, pyo.Var):
                    if (options.ncp_finite > 1
                            and hasattr(m, f"{sv_name}_interpolation_constraints")):
                        getattr(m, f"{sv_name}_interpolation_constraints").deactivate()

        m.obj = pyo.Objective(expr=1)

        if plant_options.model_output_dir:
            import os
            os.makedirs(plant_options.model_output_dir, exist_ok=True)
            with open(os.path.join(plant_options.model_output_dir, "plant_model.txt"), "w") as f:
                m.pprint(ostream=f)

        with open("plant_model.txt", "w") as f:
            m.pprint(ostream=f)

        self._model = m
        self._solver = _ipopt_solver()

    def solve(self):
        """Solve the plant model in place."""
        results = self._solver.solve(self._model, tee=self.options.tee_flag)
        try:
            _check_optimal(results, "plant")
        except RuntimeError:
            if self.options.debug_flag:
                _report_constraint_violations(self._model, label="plant failure")
            raise

    def __getattr__(self, name: str):
        # Delegate unknown attribute lookups to the underlying Pyomo model.
        # __getattr__ is only called when normal lookup fails, so 'options',
        # '_model', '_solver', and 'solve' are all found first.
        return getattr(self._model, name)
