"""
Controller classes for finite- and infinite-horizon NMPC.
"""
import resource
import pyomo.environ as pyo
from .make_model import _make_infinite_horizon_model, _make_finite_horizon_model, _ipopt_solver, _check_optimal
from .model_equations import _get_model
from .tools.indexing_tools import _add_time_indexed_expression
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

    def solve(self):
        """Solve the controller optimisation problem in place."""
        before = resource.getrusage(resource.RUSAGE_CHILDREN)
        results = self._solver.solve(self._model, tee=self.options.tee_flag)
        _check_optimal(results)
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

        if options.custom_objective:
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
                return obj

            m.lyapunov = pyo.Expression(
                expr=m.infinite_block.phi_track[m.infinite_block.time.last()]
                * options.beta / options.sampling_time
            )

        m.objective = pyo.Objective(rule=objective_rule, sense=pyo.minimize)

        self._model = m

        with open("model_output.txt", "w") as f:
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

        if options.custom_objective:
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
                return stage_cost

        m.objective = pyo.Objective(rule=objective_rule, sense=pyo.minimize)

        self._model = m

        with open("model_output.txt", "w") as f:
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
