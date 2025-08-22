from pyomo.environ import Objective, SolverFactory, Suffix
from workspace.models.Diehl_dist.diehl_dist_model_w_t0_con import (
    dist_daemodel,
    create_model_hard_bounds,
)
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
)

def solve_dist_ss_problem(
        economic_bound_function, economic_stage_cost_function, get_duals=False,
        adjustFeed=None,
    ):
    ssMod = dist_daemodel(dynamic=False)
    create_model_hard_bounds(ssMod)
    economic_bound_function(ssMod)
    ssMod.econ_stage_cost = economic_stage_cost_function(ssMod, ssMod.t, dynamic=False)
    ssMod.obj = Objective(expr = ssMod.econ_stage_cost[ssMod.t.first()])

    if adjustFeed is not None:
        if type(adjustFeed) is float or type(adjustFeed) is int:
            ssMod.feed[21].set_value(adjustFeed)


    if get_duals:
        ### Declare suffixes
        # Ipopt bound multipliers (obtained from solution)
        ssMod.ipopt_zL_out = Suffix(direction=Suffix.IMPORT)
        ssMod.ipopt_zU_out = Suffix(direction=Suffix.IMPORT)
        # Ipopt bound multipliers (sent to solver)
        ssMod.ipopt_zL_in = Suffix(direction=Suffix.EXPORT)
        ssMod.ipopt_zU_in = Suffix(direction=Suffix.EXPORT)
        # Obtain dual solutions from first solve and send to warm start
        ssMod.dual = Suffix(direction=Suffix.IMPORT_EXPORT)

    ssdof = degrees_of_freedom(ssMod)
    ss_solver = SolverFactory("ipopt")
    ss_solver.options['halt_on_ampl_error'] = 'yes'
    ss_solver.solve(ssMod, tee=True)

    return ssMod
