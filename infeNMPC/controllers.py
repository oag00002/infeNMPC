"""
Controller classes for finite- and infinite-horizon NMPC.
"""
import resource
import pyomo.environ as pyo
from .make_model import _make_infinite_horizon_model, _make_finite_horizon_model, _ipopt_solver, _ipopt_warm_solver, _ipopt_revive_solver, _check_optimal
from .model_equations import _get_model
from .tools.indexing_tools import _add_time_indexed_expression
from .tools.debug_tools import _report_constraint_violations
from .infNMPC_options import Options


class Controller:
    """
    Base class for NMPC controllers.

    Wraps a Pyomo ``ConcreteModel`` and provides a ``solve()`` method.
    Attribute access falls through to the underlying model so that callers
    can use ``controller.interface``, ``controller.time``, etc. unchanged.

    Parameters
    ----------
    options : Options
        Simulation configuration.
    """

    def __init__(self, options: Options):
        self.options = options
        self._model = None
        self._solver = _ipopt_solver()
        self._warm_solver = _ipopt_warm_solver()
        self._revive_solver = _ipopt_revive_solver()
        self._initialized = False  # False until the first (cold) solve completes

    def solve(self):
        """Solve the controller optimization problem in place.

        The first call uses the cold-start solver (``_ipopt_solver``).  Every
        subsequent call — after a shift-and-load has placed the previous
        solution in the model as a warm start — uses ``_ipopt_warm_solver``
        with smaller ``mu_init`` and ``bound_push``/``bound_frac`` values.
        """
        solver = self._warm_solver if self._initialized else self._solver
        before = resource.getrusage(resource.RUSAGE_CHILDREN)
        results = solver.solve(self._model, tee=self.options.tee_flag)
        try:
            _check_optimal(results)
        except RuntimeError:
            self.last_tc = results.solver.termination_condition
            if self.options.debug_flag:
                label = "controller initial solve failure" if not self._initialized else "controller solve failure"
                _report_constraint_violations(self._model, label=label)
            raise
        after = resource.getrusage(resource.RUSAGE_CHILDREN)
        self.last_solve_time = (after.ru_utime - before.ru_utime) + (after.ru_stime - before.ru_stime)
        if self.options.debug_flag:
            label = "controller initial solve" if not self._initialized else "controller solve"
            _report_constraint_violations(self._model, label=label)
        self._initialized = True

    def revive_solve(self):
        """Retry a failed solve with the revive solver from the current model state.

        Uses ``_ipopt_revive_solver`` (``bound_relax_factor=1e-8``, ``tol=1e-6``),
        starting from whatever state the model is in after the failed warm solve
        (typically the restoration point, which already has tiny ``inf_pr``).
        Does not change ``_initialized`` — the warm solver remains active on the
        next MPC iteration.  ``last_solve_time`` is set to the revive solve time
        only; the failed initial solve is not counted.
        """
        before = resource.getrusage(resource.RUSAGE_CHILDREN)
        results = self._revive_solver.solve(self._model, tee=self.options.tee_flag)
        try:
            _check_optimal(results)
        except RuntimeError:
            self.last_tc = results.solver.termination_condition
            raise
        after = resource.getrusage(resource.RUSAGE_CHILDREN)
        self.last_solve_time = (after.ru_utime - before.ru_utime) + (after.ru_stime - before.ru_stime)

    def __getattr__(self, name: str):
        return getattr(self._model, name)


class InfiniteHorizonController(Controller):
    """
    Infinite-horizon NMPC controller.

    Builds an infinite-horizon Pyomo model, attaches the objective, performs
    an initial solve, and writes the model structure to ``model_output.txt``.

    Parameters
    ----------
    options : Options
        Controller configuration including weights, solver options, and
        custom objective flag.
    data : optional
        Initial data loaded into the model before the first solve.
    """

    def __init__(self, options: Options, data=None):
        super().__init__(options)

        m = pyo.ConcreteModel()
        m = _make_infinite_horizon_model(m, options)

        if data is not None:
            m.interface.load_data(data, time_points=list(m.time))

        print('Generating Infinite Horizon Controller')

        m.finite_block.stage_cost_index = pyo.Set(
            initialize=lambda m: list(m.CV_index) + list(m.MV_index)
        )

        # ---- Objective ----
        _use_soft_tc = options.terminal_constraint_type == 'soft'
        _use_soft_lyap = (options.lyap_flag
                          and getattr(options, 'lyap_constraint_type', 'hard') == 'soft')

        if options.objective == 'economic':
            custom_objective = _get_model(options).custom_objective
            cost_fn = custom_objective(m.finite_block, options)

            def objective_rule(m):
                finite_elements = m.finite_block.time.get_finite_elements()
                stage_cost = sum(
                    cost_fn(m.finite_block, t) - m.finite_block.ss_obj_value
                    for t in finite_elements
                    if t != m.finite_block.time.first()
                )
                terminal_cost = (
                    m.infinite_block.phi[m.infinite_block.time.last()]
                    * options.beta / options.sampling_time
                )
                obj = stage_cost + terminal_cost
                if _use_soft_tc:
                    obj = obj + (m.infinite_block.infinite_terminal_soft_penalty
                                 * options.terminal_soft_weight)
                if _use_soft_lyap:
                    obj = obj + m.lyap_slack * options.lyap_soft_weight
                return obj

            m.lyapunov = pyo.Expression(
                expr=m.infinite_block.phi_track[m.infinite_block.time.last()]
                * options.beta / options.sampling_time
            )

        elif options.terminal_cost_riemann:
            def objective_rule(m):
                c = options.stage_cost_weights
                finite_elements = m.finite_block.time.get_finite_elements()
                stage_cost = sum(
                    sum(
                        c[i] * (
                            _add_time_indexed_expression(m.finite_block, var_name, t)
                            - m.finite_block.steady_state_values[var_name]
                        ) ** 2
                        for i, var_name in enumerate(m.finite_block.stage_cost_index)
                    )
                    for t in finite_elements
                    if t != m.finite_block.time.first()
                )

                tau_list = list(m.infinite_block.time.data())
                tau_prev = 0
                obj_infinite = 0
                for tau in tau_list:
                    if m.infinite_block.time.first() < tau < m.infinite_block.time.last():
                        tracking = sum(
                            c[i] * (
                                _add_time_indexed_expression(m.infinite_block, var_name, tau)
                                - m.infinite_block.steady_state_values[var_name]
                            ) ** 2
                            for i, var_name in enumerate(m.infinite_block.stage_cost_index)
                        )
                        obj_infinite += tracking * (tau - tau_prev) / (
                            m.infinite_block.gamma / options.sampling_time * (1 - tau ** 2)
                        )
                        tau_prev = tau

                obj = stage_cost + obj_infinite * options.beta / options.sampling_time
                if _use_soft_tc:
                    obj = obj + (m.infinite_block.infinite_terminal_soft_penalty
                                 * options.terminal_soft_weight)
                if _use_soft_lyap:
                    obj = obj + m.lyap_slack * options.lyap_soft_weight
                if hasattr(m.finite_block, 'slack_index') and options.slack_penalty_weight > 0:
                    obj += options.slack_penalty_weight * sum(
                        getattr(m.finite_block, sv_name)[idx]
                        for sv_name in m.finite_block.slack_index
                        for idx in getattr(m.finite_block, sv_name).index_set()
                    )
                if hasattr(m.infinite_block, 'slack_index') and options.slack_penalty_weight > 0:
                    obj += options.slack_penalty_weight * sum(
                        getattr(m.infinite_block, sv_name)[idx]
                        for sv_name in m.infinite_block.slack_index
                        for idx in getattr(m.infinite_block, sv_name).index_set()
                    )
                return obj

        else:
            def objective_rule(m):
                c = options.stage_cost_weights
                finite_elements = m.finite_block.time.get_finite_elements()
                stage_cost = sum(
                    sum(
                        c[i] * (
                            _add_time_indexed_expression(m.finite_block, var_name, t)
                            - m.finite_block.steady_state_values[var_name]
                        ) ** 2
                        for i, var_name in enumerate(m.finite_block.stage_cost_index)
                    )
                    for t in finite_elements
                    if t != m.finite_block.time.first()
                )

                if options.input_suppression:
                    for t in m.finite_block.time:
                        if t > m.finite_block.time.first():
                            t_prev = m.finite_block.time.prev(t)
                            for mv in m.finite_block.MV_index:
                                mv_expr = (
                                    _add_time_indexed_expression(m.finite_block, mv, t)
                                    - _add_time_indexed_expression(m.finite_block, mv, t_prev)
                                )
                                stage_cost += options.input_suppression_factor * mv_expr ** 2

                    for t in m.infinite_block.time:
                        if t > m.infinite_block.time.first():
                            t_prev = m.infinite_block.time.prev(t)
                            for mv in m.infinite_block.MV_index:
                                mv_expr = (
                                    _add_time_indexed_expression(m.infinite_block, mv, t)
                                    - _add_time_indexed_expression(m.infinite_block, mv, t_prev)
                                )
                                stage_cost += options.input_suppression_factor * mv_expr ** 2

                terminal_cost = (
                    m.infinite_block.phi[m.infinite_block.time.last()]
                    * options.beta / options.sampling_time
                )
                obj = stage_cost + terminal_cost
                if _use_soft_tc:
                    obj = obj + (m.infinite_block.infinite_terminal_soft_penalty
                                 * options.terminal_soft_weight)
                if _use_soft_lyap:
                    obj = obj + m.lyap_slack * options.lyap_soft_weight
                if hasattr(m.finite_block, 'slack_index') and options.slack_penalty_weight > 0:
                    obj += options.slack_penalty_weight * sum(
                        getattr(m.finite_block, sv_name)[idx]
                        for sv_name in m.finite_block.slack_index
                        for idx in getattr(m.finite_block, sv_name).index_set()
                    )
                if hasattr(m.infinite_block, 'slack_index') and options.slack_penalty_weight > 0:
                    obj += options.slack_penalty_weight * sum(
                        getattr(m.infinite_block, sv_name)[idx]
                        for sv_name in m.infinite_block.slack_index
                        for idx in getattr(m.infinite_block, sv_name).index_set()
                    )
                return obj

            m.lyapunov = pyo.Expression(
                expr=m.infinite_block.phi_track[m.infinite_block.time.last()]
                * options.beta / options.sampling_time
            )

        m.objective = pyo.Objective(rule=objective_rule, sense=pyo.minimize)

        self._model = m

        if options.model_output_dir:
            import os
            os.makedirs(options.model_output_dir, exist_ok=True)
            out_path = os.path.join(options.model_output_dir, "controller_model.txt")
        else:
            out_path = "model_output.txt"
        with open(out_path, "w") as f:
            m.pprint(ostream=f)

        print('Infinite Horizon Controller Initial Solve')
        self.solve()


class FiniteHorizonController(Controller):
    """
    Finite-horizon NMPC controller.

    Builds a finite-horizon Pyomo model, attaches the objective, performs
    an initial solve, and writes the model structure to ``model_output.txt``.

    Parameters
    ----------
    options : Options
        Controller configuration including weights, solver options, and
        custom objective flag.
    data : optional
        List of initial data objects loaded into the model before the first solve.
    """

    def __init__(self, options: Options, data=None):
        super().__init__(options)

        m = pyo.ConcreteModel()
        m = _make_finite_horizon_model(m, options)

        if data is not None:
            print("Loading initial data into finite horizon controller")
            for i in range(len(data)):
                m.interface.load_data(data[i], time_points=list(m.time)[i])

        print('Generating Finite Horizon Controller')

        m.stage_cost_index = pyo.Set(
            initialize=lambda m: list(m.CV_index) + list(m.MV_index)
        )

        # ---- Objective ----
        _use_soft_tc = options.terminal_constraint_type == 'soft'
        _use_soft_lyap = (options.lyap_flag
                          and getattr(options, 'lyap_constraint_type', 'hard') == 'soft')

        if options.objective == 'economic':
            custom_objective = _get_model(options).custom_objective
            cost_fn = custom_objective(m, options)

            def objective_rule(m):
                finite_elements = m.time.get_finite_elements()
                stage_cost = sum(
                    cost_fn(m, t) - m.ss_obj_value
                    for t in finite_elements
                    if t != m.time.first()
                )
                if _use_soft_tc:
                    stage_cost = (stage_cost
                                  + m.finite_terminal_soft_penalty
                                  * options.terminal_soft_weight)
                if _use_soft_lyap:
                    stage_cost = stage_cost + m.lyap_slack * options.lyap_soft_weight
                return stage_cost

        else:
            def objective_rule(m):
                c = options.stage_cost_weights
                finite_elements = m.time.get_finite_elements()
                stage_cost = sum(
                    sum(
                        c[i] * (
                            _add_time_indexed_expression(m, var_name, t)
                            - m.steady_state_values[var_name]
                        ) ** 2
                        for i, var_name in enumerate(m.stage_cost_index)
                    )
                    for t in finite_elements
                    if t != m.time.first()
                )

                if options.input_suppression:
                    for t in m.time:
                        if t > m.time.first():
                            t_prev = m.time.prev(t)
                            for mv in m.MV_index:
                                mv_expr = (
                                    _add_time_indexed_expression(m, mv, t)
                                    - _add_time_indexed_expression(m, mv, t_prev)
                                )
                                stage_cost += options.input_suppression_factor * mv_expr ** 2

                if _use_soft_tc:
                    stage_cost = (stage_cost
                                  + m.finite_terminal_soft_penalty
                                  * options.terminal_soft_weight)
                if _use_soft_lyap:
                    stage_cost = stage_cost + m.lyap_slack * options.lyap_soft_weight
                if hasattr(m, 'slack_index') and options.slack_penalty_weight > 0:
                    stage_cost += options.slack_penalty_weight * sum(
                        getattr(m, sv_name)[idx]
                        for sv_name in m.slack_index
                        for idx in getattr(m, sv_name).index_set()
                    )
                return stage_cost

        m.objective = pyo.Objective(rule=objective_rule, sense=pyo.minimize)

        self._model = m

        if options.model_output_dir:
            import os
            os.makedirs(options.model_output_dir, exist_ok=True)
            out_path = os.path.join(options.model_output_dir, "controller_model.txt")
        else:
            out_path = "model_output.txt"
        with open(out_path, "w") as f:
            m.pprint(ostream=f)

        print('Finite Horizon Controller Initial Solve')
        self.solve()


# ---------------------------------------------------------------------------
# Backwards-compatible factory aliases
# ---------------------------------------------------------------------------

def _make_infinite_horizon_controller(options: Options, data=None) -> InfiniteHorizonController:
    """Backwards-compatible factory wrapper for ``InfiniteHorizonController``."""
    return InfiniteHorizonController(options, data=data)


def _make_finite_horizon_controller(options: Options, data=None) -> FiniteHorizonController:
    """Backwards-compatible factory wrapper for ``FiniteHorizonController``."""
    return FiniteHorizonController(options, data=data)
