"""
Plant model for closed-loop NMPC simulations.
"""
import pyomo.environ as pyo
from .make_model import _make_finite_horizon_model, _ipopt_warm_solver, _check_optimal
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

        # Spec slack vars (spec_slack_index) — fixed at 0 in the plant.
        # The plant is a pure physics simulation (obj=1 constant); spec constraints
        # have no penalty so spec slacks are genuinely free, making the NLP
        # under-determined (9 extra DOF for 3 slacks × 3 collocation points).
        # Fix them at 0 and deactivate the associated spec lb constraints so the
        # plant NLP is square.  Models advertise the constraint names via the
        # optional spec_lb_constraint_index Set.
        if hasattr(m, 'spec_slack_index'):
            for sv_name in m.spec_slack_index:
                sv_var = getattr(m, sv_name, None)
                if sv_var is not None and isinstance(sv_var, pyo.Var):
                    if (options.ncp_finite > 1
                            and hasattr(m, f"{sv_name}_interpolation_constraints")):
                        getattr(m, f"{sv_name}_interpolation_constraints").deactivate()
                    sv_var.fix(0)
        if hasattr(m, 'spec_lb_constraint_index'):
            for con_name in m.spec_lb_constraint_index:
                con = getattr(m, con_name, None)
                if con is not None:
                    con.deactivate()

        m.obj = pyo.Objective(expr=1)

        if plant_options.model_output_dir:
            import os
            os.makedirs(plant_options.model_output_dir, exist_ok=True)
            with open(os.path.join(plant_options.model_output_dir, "plant_model.txt"), "w") as f:
                m.pprint(ostream=f)

        if plant_options.debug_flag:
            with open("plant_model.txt", "w") as f:
                m.pprint(ostream=f)

        self._model = m
        # Plant is always warm-started from the controller's first-FE solution,
        # so it never needs strict bound enforcement.  _ipopt_warm_solver includes
        # bound_relax_factor=1e-8 and acceptable_* settings that prevent the same
        # "restoration phase converged to small primal infeasibility" pattern seen
        # in the controller warm solves.
        self._solver = _ipopt_warm_solver()

    def solve(self):
        """Solve the plant model in place."""
        results = self._solver.solve(self._model, tee=self.options.tee_flag)
        try:
            _check_optimal(results, "plant")
        except RuntimeError:
            if self.options.debug_flag:
                _report_constraint_violations(self._model, label="plant failure")
            raise
        if self.options.debug_flag:
            _report_constraint_violations(self._model, label="plant solve")

    def __getattr__(self, name: str):
        # Delegate unknown attribute lookups to the underlying Pyomo model.
        # __getattr__ is only called when normal lookup fails, so 'options',
        # '_model', '_solver', and 'solve' are all found first.
        return getattr(self._model, name)
