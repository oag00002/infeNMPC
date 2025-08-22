from pyomo.environ import (
    Set,
    Var,
    Param,
    Expression,
    Constraint,
    SolverFactory,
#     Objective,
    ComponentUID,
    Reference,
    Block
)
from pyomo.environ import value as pyo_val

from workspace.models.Diehl_dist.diehl_dist_model_w_t0_con import (
    dist_daemodel,
#     create_model_hard_bounds,
)
from workspace.models.Diehl_dist.diehl_dist_econ_stage_cost_bound import (
    dist_econ_stage_cost,
    dist_econ_bound,
)
from workspace.models.Diehl_dist.dist_ss_problem import (
    solve_dist_ss_problem
)
from workspace.dynamics import (
    discretize_and_create_fe_set,
    initialize_vars,
    retrieve_dae_vars_at_tp
)
from workspace.solve_NLP import(
    solve_problem,
)
from workspace.eNMPC_red_grad.Diehl_dist.\
    self_stabilizing_KKT_expression_alternate_stage_cost import(
        construct_KKT_expressions,
        initialize_KKT_variables_with_ipopt,
        construct_KKT_constraints_for_stages,
        initialize_dual_variables_by_file,
)
import numpy as np
import time
import json
from os.path import join
# import matplotlib.pyplot as plt

initialize_KKT_variable_method = 2 #or2
economic_stage_cost_num = 1
economic_stage_cost_fun = dist_econ_stage_cost(economic_stage_cost_num)
economic_bound_fun = dist_econ_bound(economic_stage_cost_num)

str_states = ("x", "M")
str_controls = ("Rec", "Qr")

# ipopt solver
csolver = SolverFactory("ipopt")
csolver.options['linear_solver'] = 'ma57'
csolver.options['bound_push'] = 1.0E-8

ssMod = solve_dist_ss_problem(economic_bound_fun, economic_stage_cost_fun)
ss_data = retrieve_dae_vars_at_tp(ssMod, ssMod.t.first())

# time0 = time.time()

N = 10
h = 60. #sampling time
olnmpc = dist_daemodel(nfe=N, th=h)
olnmpc.name = "self_stabilizing enmpc_dist"

discretize_and_create_fe_set(olnmpc,
                             nfe=N,
                             ncp=3,
                             is_control_model=True,
                             control_vars=[olnmpc.Rec, olnmpc.Qr],
                             )

initialize_vars(olnmpc, ss_data)


# Economic stage cost expression
olnmpc.econ_stage_cost = economic_stage_cost_fun(olnmpc, olnmpc.Nfe, dynamic=True)


if initialize_KKT_variable_method == 1:
    construct_KKT_expressions(olnmpc)
    initialize_KKT_variables_with_ipopt(olnmpc, csolver)
    construct_KKT_constraints_for_stages(olnmpc)

elif initialize_KKT_variable_method == 2:

    construct_KKT_expressions(olnmpc)
    construct_KKT_constraints_for_stages(olnmpc)
    initialize_dual_variables_by_file(olnmpc)

    for var in olnmpc.component_objects(Var):
        try:
            var.fix()
        except AttributeError:
            var[...].fix()
    for con in olnmpc.component_objects(Constraint):
        con.deactivate()

    olnmpc.KKT1_var.unfix()
    olnmpc.KKT2_var.unfix()
    olnmpc.KKT3_var.unfix()
    olnmpc.KKT4_var.unfix()

    olnmpc.KKT1.activate()
    olnmpc.KKT2.activate()
    olnmpc.KKT3.activate()
    olnmpc.KKT4.activate()

    solve_problem("initialize KKT vars", csolver, olnmpc)

    for var in olnmpc.component_objects(Var):
        try:
            var.unfix()
        except AttributeError:
            var[...].unfix()
    for con in olnmpc.component_objects(Constraint):
        con.activate()
