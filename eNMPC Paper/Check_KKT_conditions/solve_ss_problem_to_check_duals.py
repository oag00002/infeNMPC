from workspace.models.Diehl_dist.dist_ss_problem import (
    solve_dist_ss_problem
)
from workspace.models.Diehl_dist.diehl_dist_econ_stage_cost_bound import (
    dist_econ_stage_cost,
    dist_econ_bound,
)

economic_stage_cost_num = 1
economic_stage_cost_fun = dist_econ_stage_cost(economic_stage_cost_num)
economic_bound_fun = dist_econ_bound(economic_stage_cost_num)

ssMod = solve_dist_ss_problem(
    economic_bound_fun, economic_stage_cost_fun, get_duals=True)

for i, j in ssMod.ipopt_zL_out.items():
    if abs(j) >= 1.0E-4:
        print(i, j)

for i, j in ssMod.ipopt_zU_out.items():
    if abs(j) >= 1.0E-4:
        print(i, j)
