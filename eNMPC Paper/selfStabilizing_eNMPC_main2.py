from pyomo.environ import (
    Var,
    Param,
    Constraint,
    SolverFactory,
    Objective,
    Expression,
    ComponentUID,
    Reference,
    TerminationCondition,
)
from pyomo.environ import value as pyo_val
from pyomo.dae.flatten import flatten_dae_components
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
)
from workspace.models.Diehl_dist.diehl_dist_model_w_t0_con import (
    dist_daemodel,
    create_model_hard_bounds,
)
from workspace.models.Diehl_dist.dist_ss_problem import (
    solve_dist_ss_problem
)
from workspace.models.Diehl_dist.diehl_dist_econ_stage_cost_bound import (
    dist_econ_stage_cost,
    dist_econ_bound,
    dist_tighter_control_bound_for_eNMPC,
)
from workspace.eNMPC_red_grad.Diehl_dist.\
    self_stabilizing_KKT_expression_alternate_stage_cost import(
        construct_KKT_expressions,
        initialize_KKT_variables_with_ipopt,
        construct_KKT_constraints_for_stages,
        initialize_dual_variables_by_file,
)
from workspace.dynamics import (
    discretize_and_create_fe_set,
    shift_dae_vars_by_dt,
    construct_L1relaxed_endpoint_cons,
    construct_L2relaxed_endpoint_cons,
    set_initial_conds,
    initialize_vars,
    get_comp_value_at_tp,
)
from workspace.save_dynamic_data import (
    initialize_record_comp_values,
    record_comp_values,
    save_dict_in_json,
)
from workspace.solve_NLP import(
    solve_problem,
)

import json
from os.path import join, exists
from os import makedirs
import matplotlib.pyplot as plt

import time


# Options
iterations = 45 #100
delta_sc = 1.0 #0.8 #1.0 #2.0
economic_stage_cost_num = 1
has_endpoint_constraint = True
if has_endpoint_constraint:
    endpoint_constraint_type = "L2"
endpt_con_is_active_at_step1 = True
initial_state_dir = \
    "/home/khl/KHL/workspace/workspace/models/Diehl_dist/ic_and_ss"
    # "/home/kuanhanl/KHL/Research/workspace/workspace/models/Diehl_dist/ic_and_ss"
# initial_state_filename = "dist_complicated_ic.json"
# initial_state_filename = "100iters_after_complicated_ic.json"
initial_state_filename = "13iters_after_complicated_ic.json"
initial_state_path = join(initial_state_dir, initial_state_filename)

save_result_dir = \
    "/home/khl/KHL/Distillation_column/eNMPC_sc_gradLagran/Case6_L2EndpointPenalty_rerun"

save_each_sol = True
if save_each_sol:
    save_each_sol_folder = join(save_result_dir, "sol_each_iter")
    if not exists(save_each_sol_folder):
        print("Create a folder to save the solution at each iteration!")
        makedirs(save_each_sol_folder)

economic_stage_cost_fun = dist_econ_stage_cost(economic_stage_cost_num)
economic_bound_fun = dist_econ_bound(economic_stage_cost_num)
economic_tighter_control_bound_fun = \
    dist_tighter_control_bound_for_eNMPC(economic_stage_cost_num)

str_states = ("x", "M")
str_controls = ("Rec", "Qr")


# Some functions --------------------------------------------------------------


# Load initial and steady states for setpoints --------------------------------
with open(initial_state_path, "r") as fic:
    initial_states = json.load(fic)

ssMod = solve_dist_ss_problem(economic_bound_fun, economic_stage_cost_fun)

ss_info = dict()
ss_t0 = ssMod.t.first()
for str_state in str_states:
    state = ssMod.find_component(str_state)
    for j in ssMod.tray:
        ss_info[ComponentUID(state[:, j])] = state[ss_t0, j].value
for str_control in str_controls:
    control = ssMod.find_component(str_control)
    ss_info[ComponentUID(control)] = control[ss_t0].value


# Ipopt solver ---------------------------------------------------------------
csolver = SolverFactory("ipopt")
csolver.options['linear_solver'] = 'ma57'
csolver.options['bound_push'] = 1.0E-8
# csolver.options['max_iter'] = 500

psolver = SolverFactory("ipopt")
psolver.options['linear_solver'] = 'ma57'
# psolver.options['halt_on_ampl_error'] = 'yes'


# Create NMPC ----------------------------------------------------------------
N = 18 #19 #20 #10 #15
h = 60. #sampling time
olnmpc = dist_daemodel(nfe=N, th=h)
olnmpc.name = "economic nmpc_dist"

discretize_and_create_fe_set(olnmpc,
                             nfe=N,
                             ncp=3,
                             is_control_model=True,
                             control_vars=[olnmpc.Rec, olnmpc.Qr],
                             )

cflatten_states = [Reference(
            olnmpc.find_component(ComponentUID(str_state))[:, j]
            ) for str_state in str_states for j in olnmpc.tray]
cflatten_controls = [olnmpc.find_component(ComponentUID(str_control))
                      for str_control in str_controls]

create_model_hard_bounds(olnmpc)

# Use L2 penalty on relaxed bounds!
olnmpc.x42_lbVar, olnmpc.x42_lbCon, olnmpc.x1_ubVar, olnmpc.x1_ubCon = \
    economic_bound_fun(olnmpc, dynamic=True, relaxBounds=True)

olnmpc.Rec_ubVar, olnmpc.Rec_ubCon, olnmpc.Qr_ubVar, olnmpc.Qr_ubCon = \
    economic_tighter_control_bound_fun(olnmpc)
# olnmpc.D[:].setlb(17.8) Should I use this?? No!

# NMPC economic objective
olnmpc.econ_stage_cost = economic_stage_cost_fun(olnmpc, olnmpc.Nfe, dynamic=True)

# Construct KKT expressions
construct_KKT_expressions(olnmpc)
# initialize_KKT_variables_with_ipopt(olnmpc, csolver)
construct_KKT_constraints_for_stages(olnmpc)
initialize_dual_variables_by_file(olnmpc)

# Construct relaxed stabilizing constraint
olnmpc.pre_LyaFun = Param(initialize=1.0E3, mutable=True)
olnmpc.pre_stageCost = Param(initialize=1.0E2, mutable=True)
olnmpc.cur_LyaFun = Var(initialize=1.0E3)
olnmpc.sc_relaxVar = Var(initialize=0.0, bounds=(0.0, None))
olnmpc.delta_sc = Param(initialize=delta_sc, mutable=True)
Qmat = {str_states[0]: 1.0E1,  #1.0E3
        str_states[1]: 1.0E-2} #1.0
KKTmat = {"KKT1":1.0E1,
          "KKT2":1.0E-2,
          "KKT3":1.0E0,
          "KKT4":1.0E0}
def stabilizing_con_KKT_alternate_stageCost_rule(m):
    return m.cur_LyaFun - m.pre_LyaFun <= \
            -olnmpc.delta_sc*m.pre_stageCost + m.sc_relaxVar
olnmpc.stabilizing_con_KKT_alternate_stageCost = Constraint(
    rule=stabilizing_con_KKT_alternate_stageCost_rule)
def cur_LyaFun_con_rule(m):
    return m.cur_LyaFun == sum(
        sum(
            sum(
                Qmat[str_state]*(
                    olnmpc.find_component(str_state)[i+olnmpc.h, j]-\
                        olnmpc.find_component(str_state)[i, j]
                        )**2 for j in olnmpc.tray)
            for str_state in str_states) + \
        sum(KKTmat["KKT1"]*olnmpc.KKT1_var[i, k]**2 for k in olnmpc.x_ind) + \
        sum(KKTmat["KKT2"]*olnmpc.KKT2_var[i, k]**2 for k in olnmpc.y_ind) + \
        sum(KKTmat["KKT3"]*olnmpc.KKT3_var[i, k]**2 for k in olnmpc.u_ind) + \
        sum(KKTmat["KKT4"]*olnmpc.KKT4_var[i, k]**2 for k in olnmpc.h_ind)
                                for i in olnmpc.Nfe if i!=olnmpc.Nfe.last())
olnmpc.cur_LyaFun_con = Constraint(rule=cur_LyaFun_con_rule)

if has_endpoint_constraint:
    steadystates = [ss_info[ComponentUID(state.referent)] for state in cflatten_states]
    if endpoint_constraint_type == "L1":
        raise RuntimeError("L1 penalty is not supported!")
        # olnmpc.rho1 = Param(initialize=1.0E8, mutable=True)
        # olnmpc.rho2 = Param(initialize=1.0E4, mutable=True)
        # olnmpc.endpt_set, olnmpc.endpt_varTypeSet, olnmpc.endpt_var, olnmpc.endpt_con = \
        #     construct_L1relaxed_endpoint_cons(cflatten_states, steadystates, olnmpc.t.last())


        # olnmpc.econ_obj = Objective(expr = sum(olnmpc.econ_stage_cost[i]
        #                                            for i in olnmpc.Nfe if i!=olnmpc.t.last())+
        #                                         olnmpc.rho1*sum(olnmpc.endpt_var[k,"p"] + olnmpc.endpt_var[k,"n"]
        #                                                         for k in range(0, 42))+
        #                                         olnmpc.rho2*sum(olnmpc.endpt_var[k,"p"] + olnmpc.endpt_var[k,"n"]
        #                                                         for k in range(42, 84))+
        #                                         olnmpc.rho1*sum(olnmpc.x42_lbVar[i]**2 + olnmpc.x1_ubVar[i]**2
        #                                                         for i in olnmpc.t if i>0)+
        #                                         olnmpc.rho2*sum(olnmpc.Rec_ubVar[i]**2 + olnmpc.Qr_ubVar[i]**2
        #                                                         for i in olnmpc.t if i>0)
        #                                             )
    elif endpoint_constraint_type == "L2":
        olnmpc.rho1 = Param(initialize=1.0E4, mutable=True) #1.0E8 1.0E4
        olnmpc.rho2 = Param(initialize=1.0E2, mutable=True) #1.0E4 1.0E2
        olnmpc.rho3 = Param(initialize=1.0E4, mutable=True)

        olnmpc.rho4 = Param(initialize=1.0E4, mutable=True)
        olnmpc.rho5 = Param(initialize=1.0E2, mutable=True)

        olnmpc.state_endptVar = Var(olnmpc.x_ind, initialize=0.0)
        def state_endpt_con_rule(m, i):
            tf = m.Nfe.last()
            tf2 = m.Nfe.at(-2)
            return cflatten_states[i][tf] - cflatten_states[i][tf2] == m.state_endptVar[i]
        olnmpc.state_endpt_con = Constraint(olnmpc.x_ind, rule=state_endpt_con_rule)

        tf2 = olnmpc.Nfe.at(-2)
        olnmpc.econ_obj = Objective(expr = sum(olnmpc.econ_stage_cost[i]
                                                   for i in olnmpc.Nfe if i!=olnmpc.t.last()) +
                                                olnmpc.rho3*olnmpc.sc_relaxVar**2 +
                                                olnmpc.rho4*sum(olnmpc.state_endptVar[k]**2
                                                                for k in range(0, 42)) +
                                                olnmpc.rho5*sum(olnmpc.state_endptVar[k]**2
                                                                for k in range(42, 84)) +
                                                olnmpc.rho4*sum(olnmpc.KKT1_var[tf2, k]**2
                                                                for k in olnmpc.x_ind) +
                                                olnmpc.rho4*sum(olnmpc.KKT2_var[tf2, k]**2
                                                                for k in olnmpc.y_ind) +
                                                olnmpc.rho4*sum(olnmpc.KKT3_var[tf2, k]**2
                                                                for k in olnmpc.u_ind) +
                                                olnmpc.rho4*sum(olnmpc.KKT4_var[tf2, k]**2
                                                                for k in olnmpc.h_ind) +
                                                olnmpc.rho1*sum(olnmpc.x42_lbVar[i]**2 + olnmpc.x1_ubVar[i]**2
                                                                for i in olnmpc.t if i>0) +
                                                olnmpc.rho2*sum(olnmpc.Rec_ubVar[i]**2 + olnmpc.Qr_ubVar[i]**2
                                                                for i in olnmpc.t if i>0)
                                                    )
    else:
        raise RuntimeError("Declared endpoint penalty type is not supported!")
else:
    olnmpc.rho1 = Param(initialize=1.0E4, mutable=True) #1.0E8
    olnmpc.rho2 = Param(initialize=1.0E2, mutable=True) #1.0E4
    olnmpc.rho3 = Param(initialize=1.0E4, mutable=True)
    olnmpc.econ_obj = Objective(expr = sum(olnmpc.econ_stage_cost[i]
                                               for i in olnmpc.Nfe if i!=olnmpc.t.last()) +
                                            olnmpc.rho3*olnmpc.sc_relaxVar**2 +
                                            olnmpc.rho1*sum(olnmpc.x42_lbVar[i]**2 + olnmpc.x1_ubVar[i]**2
                                                for i in olnmpc.t if i>0) +
                                            olnmpc.rho2*sum(olnmpc.Rec_ubVar[i]**2 + olnmpc.Qr_ubVar[i]**2
                                                for i in olnmpc.t if i>0)
                                    )

# Initialize vars
initialize_vars(olnmpc, initial_states)

set_initial_conds(olnmpc.t.first(), cflatten_states, initial_states)
olnmpc.Rec[0].fix()
olnmpc.Qr[0].fix()

cscalar_vars, cdae_vars = flatten_dae_components(olnmpc, olnmpc.t, Var)
cscalar_cons, cdae_cons = flatten_dae_components(olnmpc, olnmpc.t, Constraint)


# Create Plant-----------------------------------------------------------------
plant = dist_daemodel(nfe=1, th=h)
plant.name = "plant"

discretize_and_create_fe_set(plant,
                             nfe=1,
                             ncp=3,
                             is_control_model=False,
                             )

pflatten_states = [Reference(
            plant.find_component(ComponentUID(str_state))[:, j]
            ) for str_state in str_states for j in plant.tray]
pflatten_controls = [plant.find_component(ComponentUID(str_control))
                      for str_control in str_controls]

plant.econ_stage_cost = economic_stage_cost_fun(plant, plant.Nfe, dynamic=True)

# Initialize vars
initialize_vars(plant, initial_states)

# plant.x[:,42].setlb(0.98)
# plant.x[0,42].setlb(0.0) # Relax the bound at t0

set_initial_conds(plant.t.first(), pflatten_states, initial_states)
plant.Rec.fix()
plant.Qr.fix()

pscalar_vars, pdae_vars = flatten_dae_components(plant, plant.t, Var)


# Start iterations -----------------------------------------------------------
states_of_interest = pflatten_states + [plant.D, Reference(plant.L[:, 1])]
state_record = initialize_record_comp_values(states_of_interest, init_t=plant.t.first())
economic_cost_record = initialize_record_comp_values([plant.econ_stage_cost], init_t=None)
control_record = initialize_record_comp_values(cflatten_controls, init_t=None)
cur_LyaFun_record = initialize_record_comp_values([olnmpc.cur_LyaFun], init_t=None)
pre_stageCost_record = initialize_record_comp_values([olnmpc.pre_stageCost], init_t=None)
notOptimal_step_record = {"notOptimal_step" : []}
computational_time = {"computational_time": []}

for step in range(1, iterations+1):
    print(step)

    # NMPC
    # Set ics
    if step == 1:
        # Deactivate stabilizing constraint
        olnmpc.stabilizing_con_KKT_alternate_stageCost.deactivate()
        olnmpc.cur_LyaFun_con.deactivate()
        olnmpc.sc_relaxVar.fix(0.0)
        olnmpc.cur_LyaFun.fix(1.0E3)

        # Deactivate and fix KKT vars and cons
        olnmpc.KKT1.deactivate()
        olnmpc.KKT2.deactivate()
        olnmpc.KKT3.deactivate()
        olnmpc.KKT4.deactivate()
        olnmpc.KKT1_var[...].fix(0.0)
        olnmpc.KKT2_var[...].fix(0.0)
        olnmpc.KKT3_var[...].fix(0.0)
        olnmpc.KKT4_var[...].fix(0.0)
        olnmpc.lamb[...].fix()
        olnmpc.nue[...].fix()
        olnmpc.eta[...].fix()

        if has_endpoint_constraint:
            # Deactivate and fix endpt cons and relaxed vars
            olnmpc.state_endpt_con.deactivate()
            olnmpc.state_endptVar[...].fix(0.0)

            if endpt_con_is_active_at_step1:
                olnmpc.state_endpt_con.activate()
                olnmpc.state_endptVar[...].unfix()

                olnmpc.KKT1[tf2, :].activate()
                olnmpc.KKT2[tf2, :].activate()
                olnmpc.KKT3[tf2, :].activate()
                olnmpc.KKT4[tf2, :].activate()
                olnmpc.KKT1_var[tf2, :].unfix()
                olnmpc.KKT2_var[tf2, :].unfix()
                olnmpc.KKT3_var[tf2, :].unfix()
                olnmpc.KKT4_var[tf2, :].unfix()

        set_initial_conds(olnmpc.t.first(), cflatten_states, initial_states)

    else:
        olnmpc.stabilizing_con_KKT_alternate_stageCost.activate()
        olnmpc.cur_LyaFun_con.activate()
        olnmpc.sc_relaxVar.unfix()
        olnmpc.cur_LyaFun.unfix()

        shift_dae_vars_by_dt(cdae_vars, olnmpc.t, dt=h)
        set_initial_conds(olnmpc.t.first(), cflatten_states, curr_states)

    # This is weird that ipopt read the correct dof(i.e., 140), but
    # degrees_of_freedom function read the wrong dof(20), and it passes the assert...
    # Anyway I comment this seciton for dof check.
    # if has_endpoint_constraint and endpoint_constraint_type == "L1":
    #     pen_val = pyo_val(
    #         olnmpc.rho1*sum(olnmpc.endpt_var[k,"p"] + olnmpc.endpt_var[k,"n"]
    #                         for k in range(0, 42))+
    #         olnmpc.rho2*sum(olnmpc.endpt_var[k,"p"] + olnmpc.endpt_var[k,"n"]
    #                         for k in range(42, 84))
    #         )
    #     if abs(pen_val) <= 1.0E-2 and step > 1:
    #         # In this case, we don't need endpoint constraints
    #         olnmpc.endpt_var.fix(0.0)
    #         olnmpc.endpt_con.deactivate()
    #         assert degrees_of_freedom(olnmpc) == 2*N
    #     else:
    #         assert degrees_of_freedom(olnmpc) == 2*N + len(cflatten_states)
    # else:
    #     assert degrees_of_freedom(olnmpc) == 2*N
    
    # if step in [15, 32, 39, 55]:
    #     csolver.options["max_iter"] = 1

    # solve_problem(step, csolver, olnmpc)
    print("")
    print("iteration = ", step, ", Solve " + olnmpc.name)
    start_time = time.time()
    olnmpcResults = csolver.solve(olnmpc, tee=True, load_solutions=False)
    end_time = time.time()
    
    # if step in [15, 32, 39, 55]:
    #     csolver.options.pop("max_iter")

    if olnmpcResults.solver.termination_condition == TerminationCondition.optimal:
        olnmpc.solutions.load_from(olnmpcResults)

        curr_controls = get_comp_value_at_tp(cflatten_controls, tp=h) #get control input
        control_record = record_comp_values(
            mod=olnmpc, tp=h, record=control_record
        )

        computational_time["computational_time"].append(end_time-start_time)
        if save_each_sol:
            sol_filename = "sol_iter" + str(step) + ".json"
            sol_info = dict()
            # save vars and objective
            for comp in olnmpc.component_data_objects(ctype=[Var, Objective]):
                key = str(ComponentUID(comp))
                val = pyo_val(comp)
                sol_info[key] = val

            # save pre_LyaFun and pre_stageCost
            for comp in [olnmpc.pre_LyaFun, olnmpc.pre_stageCost, *olnmpc.econ_stage_cost.values()]:
                key = str(ComponentUID(comp))
                val = pyo_val(comp)
                sol_info[key] = val

            save_dict_in_json(sol_info, join(save_each_sol_folder, sol_filename))

    else:
        print("Optimal solutions is not found, use the previous control actions")
        notOptimal_step_record["notOptimal_step"].append(
            [step, olnmpcResults.solver.termination_condition.name]
         )

        control_record = record_comp_values(
            mod=olnmpc, tp=t0, record=control_record
        )

        computational_time["computational_time"].append("not_optimal")

    # Calculate LyaFun and first stage cost
    if step == 1:
        for var in cdae_vars:
            var[...].fix()
        for con in cdae_cons:
            con.deactivate()
        olnmpc.Rec_interpolation_constraints.deactivate()
        olnmpc.Qr_interpolation_constraints.deactivate()

        if has_endpoint_constraint:
            olnmpc.state_endptVar[:].fix()
            olnmpc.state_endpt_con.deactivate()

        olnmpc.lamb.fix()
        olnmpc.nue.fix()
        olnmpc.eta.fix()

        olnmpc.KKT1[...].activate()
        olnmpc.KKT2[...].activate()
        olnmpc.KKT3[...].activate()
        olnmpc.KKT4[...].activate()
        olnmpc.KKT1_var[...].unfix()
        olnmpc.KKT2_var[...].unfix()
        olnmpc.KKT3_var[...].unfix()
        olnmpc.KKT4_var[...].unfix()

        solve_problem("solve KKT vars @ step 1", csolver, olnmpc)

        pre_LyaFun = pyo_val(sum(
            sum(
                sum(
                    Qmat[str_state]*(
                        olnmpc.find_component(str_state)[i+olnmpc.h, j]-\
                            olnmpc.find_component(str_state)[i, j]
                            )**2 for j in olnmpc.tray)
                for str_state in str_states) + \
            sum(KKTmat["KKT1"]*olnmpc.KKT1_var[i, k]**2 for k in olnmpc.x_ind) + \
            sum(KKTmat["KKT2"]*olnmpc.KKT2_var[i, k]**2 for k in olnmpc.y_ind) + \
            sum(KKTmat["KKT3"]*olnmpc.KKT3_var[i, k]**2 for k in olnmpc.u_ind) + \
            sum(KKTmat["KKT4"]*olnmpc.KKT4_var[i, k]**2 for k in olnmpc.h_ind)
                                    for i in olnmpc.Nfe if i!=olnmpc.Nfe.last()))
        olnmpc.pre_LyaFun.set_value(pre_LyaFun)
        olnmpc.cur_LyaFun.set_value(pre_LyaFun)

        t0 = olnmpc.t.first()
        pre_stageCost = pyo_val(
            sum(
                sum(
                    Qmat[str_state]*(
                        olnmpc.find_component(str_state)[t0+olnmpc.h, j]-\
                            olnmpc.find_component(str_state)[t0, j]
                            )**2 for j in olnmpc.tray)
                for str_state in str_states) + \
            sum(KKTmat["KKT1"]*olnmpc.KKT1_var[t0, k]**2 for k in olnmpc.x_ind) + \
            sum(KKTmat["KKT2"]*olnmpc.KKT2_var[t0, k]**2 for k in olnmpc.y_ind) + \
            sum(KKTmat["KKT3"]*olnmpc.KKT3_var[t0, k]**2 for k in olnmpc.u_ind) + \
            sum(KKTmat["KKT4"]*olnmpc.KKT4_var[t0, k]**2 for k in olnmpc.h_ind)
                                    )
        olnmpc.pre_stageCost.set_value(pre_stageCost)


        for var in cdae_vars:
            var[...].unfix()
        olnmpc.Rec[0].fix()
        olnmpc.Qr[0].fix()
        for con in cdae_cons:
            con.activate()
        olnmpc.Rec_interpolation_constraints.activate()
        olnmpc.Qr_interpolation_constraints.activate()

        if has_endpoint_constraint:
            olnmpc.state_endptVar[:].unfix()
            olnmpc.state_endpt_con.activate()

        # !!! kind of cheating not to unfix, but it's ok for now !!!
        # olnmpc.lamb.unfix()
        # olnmpc.nue.unfix()
        # olnmpc.eta.unfix()

    else:
        olnmpc.pre_LyaFun.set_value(olnmpc.cur_LyaFun.value)
        t0 = olnmpc.t.first()
        pre_stageCost = pyo_val(
            sum(
                sum(
                    Qmat[str_state]*(
                        olnmpc.find_component(str_state)[t0+olnmpc.h, j]-\
                            olnmpc.find_component(str_state)[t0, j]
                            )**2 for j in olnmpc.tray)
                for str_state in str_states) + \
            sum(KKTmat["KKT1"]*olnmpc.KKT1_var[t0, k]**2 for k in olnmpc.x_ind) + \
            sum(KKTmat["KKT2"]*olnmpc.KKT2_var[t0, k]**2 for k in olnmpc.y_ind) + \
            sum(KKTmat["KKT3"]*olnmpc.KKT3_var[t0, k]**2 for k in olnmpc.u_ind) + \
            sum(KKTmat["KKT4"]*olnmpc.KKT4_var[t0, k]**2 for k in olnmpc.h_ind)
                                    )
        olnmpc.pre_stageCost.set_value(pre_stageCost)

    cur_LyaFun_record = record_comp_values(
        mod=olnmpc, tp=None, record=cur_LyaFun_record
            )
    pre_stageCost_record = record_comp_values(
        mod=olnmpc, tp=None, record=pre_stageCost_record
            )

    # if step>20 and step<=46:
    #     delta_sc = 6.0
    #     olnmpc.delta_sc.set_value(delta_sc)
    # elif step>46:
    #     delta_sc = 2.0
    #     olnmpc.delta_sc.set_value(delta_sc)

    print("step ", step)
    print("delta_sc", olnmpc.delta_sc.value)
    print("Current Lyapunov function: ", olnmpc.cur_LyaFun.value)
    print("Next descent amount(pre_stageCost): ", delta_sc*pyo_val(olnmpc.pre_stageCost))
    print("Upper bound of next Lyapunov function: ",
          olnmpc.cur_LyaFun.value-delta_sc*pyo_val(olnmpc.pre_stageCost))

    if delta_sc*olnmpc.pre_stageCost.value >= olnmpc.pre_LyaFun.value:
        print("")
        print("Greater descent amount than the LyaFun, reset descent amount(pre_stageCost)")
        olnmpc.pre_stageCost.set_value(olnmpc.pre_LyaFun.value/delta_sc)

        print("Current Lyapunov function: ", olnmpc.cur_LyaFun.value)
        print("Next descent amount(pre_stageCost): ", delta_sc*pyo_val(olnmpc.pre_stageCost))
        print("Upper bound of next Lyapunov function: ",
              olnmpc.cur_LyaFun.value-delta_sc*pyo_val(olnmpc.pre_stageCost))


    #Plant
    if step == 1:
        set_initial_conds(plant.t.first(), pflatten_states, initial_states)
    else:
        shift_dae_vars_by_dt(pdae_vars, plant.t, dt=h)
        set_initial_conds(plant.t.first(), pflatten_states, curr_states)

    for str_control in str_controls:
        control = plant.find_component(str_control)
        for t in list(plant.t)[1:]:
            control[t].set_value(curr_controls[str(ComponentUID(str_control))])

    solve_problem(step, psolver, plant)

    curr_states = get_comp_value_at_tp(pflatten_states, tp=h)
    state_record = record_comp_values(
        mod=plant, tp=h, record=state_record)

    economic_cost_record = record_comp_values(
        mod=plant, tp=0, record=economic_cost_record)

    if step in [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120]:
        # save vars
        vars_filename = "var_initialize_iter" + str(step+1) + ".json"
        var_info = dict()
        for comp in olnmpc.component_data_objects(Var):
            key = str(ComponentUID(comp))
            val = comp.value
            var_info[key] = val
        save_dict_in_json(var_info, vars_filename)

        # save curr_states
        curr_states_filename = "curr_state_iter" + str(step+1) + ".json"
        save_dict_in_json(curr_states, curr_states_filename)

        # save pre_LyaFun and pre_stageCost
        Lya_info = {str(ComponentUID(comp)): comp.value
                    for comp in [olnmpc.pre_LyaFun, olnmpc.pre_stageCost]}
        Lya_info_filename = "Lya_info_iter" + str(step+1) + ".json"
        save_dict_in_json(Lya_info, Lya_info_filename)


        # Save results
        path = join(save_result_dir, "state_record.json")
        save_dict_in_json(state_record, path)

        path = join(save_result_dir, "control_record.json")
        save_dict_in_json(control_record, path)

        path = join(save_result_dir, "economic_cost_record.json")
        save_dict_in_json(economic_cost_record, path)

        path = join(save_result_dir, "cur_LyaFun_record.json")
        save_dict_in_json(cur_LyaFun_record, path)

        path = join(save_result_dir, "pre_stageCost_record.json")
        save_dict_in_json(pre_stageCost_record, path)

        path = join(save_result_dir, "notOptimal_step_record.json")
        save_dict_in_json(notOptimal_step_record, path)
        
        path = join(save_result_dir, "computational_time.json")
        save_dict_in_json(computational_time, path)        
