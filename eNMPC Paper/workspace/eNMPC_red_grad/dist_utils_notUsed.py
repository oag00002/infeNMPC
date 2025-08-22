from pyomo.environ import (
    Param,
    Var,
    Constraint,
    Objective,
    SolverFactory,
    ComponentUID,
    Expression
    )
from pyomo.environ import value as pyo_val
from pyomo.dae.flatten import flatten_dae_components
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    )
from workspace.dynamics import augment_de_to_ss_eq
from workspace.models.Diehl_dist.diehl_dist_model import (
    dist_daemodel,
    create_model_hard_bounds,
    )
from workspace.models.Diehl_dist.diehl_dist_econ_stage_cost_notUsed import dist_econ_stage_cost
from workspace.eNMPC_red_grad.utils import solve_problem
import copy
import json

state_strs = ["x", "M", "T"]
control_strs = ["Rec", "Qr"]


def obtain_optimal_ss():
    h = 60
    mod = dist_daemodel(nfe=1, th=h)
    augment_de_to_ss_eq(mod)
    create_model_hard_bounds(mod)

    # # eco3 bounds
    # mod.x[:, 42].setlb(0.98)
    # mod.Rec[:].setub(1.5)

    # # eco4 bounds
    # mod.x[:, 42].setlb(0.98)
    # mod.D[:].setlb(18.0)

    # eco5 bounds
    mod.x[:, 42].setlb(0.98)
    mod.x[:, 1].setub(0.02)

    dist_econ_stage_cost(mod, mod.t, dynamic=False)
    mod.econ_obj = Objective(expr=mod.econ_stage_cost_expr[h])

    dof = degrees_of_freedom(mod)
    print("degrees of freedom for s.s. problem: ", dof)
    assert dof == 2

    solver = SolverFactory("ipopt")
    solver.solve(mod, tee = True)

    # create ss info
    tf = mod.t.last()
    scalar_var, dae_var = flatten_dae_components(mod, mod.t, ctype = Var)

    ss_info = dict()
    for var in dae_var:
        cuid = ComponentUID(var.referent)
        val = var[tf].value
        ss_info[str(cuid)] = val

    ss_info["phi_ss"] = pyo_val(mod.econ_obj)

    return ss_info, mod


def fill_vector_w_pyo_var_t0(t0, var, vec):
    nrow, ncol = vec.shape
    assert ncol == 1
    assert len([ele for ele in var[t0, :]]) == nrow

    for ind, ele in enumerate(var[t0, :]):
        vec[ind][0] = ele


def construct_relaxed_endpoint_constraints(mod, w_x, w_M, w_rg):

    # Endpoint constraint for the states
    tray_set = mod.tray
    mod.end_pen_x = Var(tray_set, initialize=0.0, bounds=(0.0, None))
    mod.end_pen_M = Var(tray_set, initialize=0.0, bounds=(0.0, None))

    def x_end_rule(m, i):
        tlast = m.Nfe.at(-1)
        tlast2 = m.Nfe.at(-2)
        return m.x[(tlast2, i)] - m.x[(tlast, i)] == m.end_pen_x[i]
    mod.x_end = Constraint(mod.tray, rule=x_end_rule)

    def M_end_rule(m, i):
        tlast = m.Nfe.at(-1)
        tlast2 = m.Nfe.at(-2)
        return m.M[(tlast2, i)] - m.M[(tlast, i)] == m.end_pen_M[i]
    mod.M_end = Constraint(mod.tray, rule=M_end_rule)

    # Endpoint constraint for the grad_u_L (reduced gradient)
    mod.end_pen_grad_u_L = Var(mod.u_ind, initialize=0.0, bounds=(0.0, None))
    def grad_u_L_end_rule(m, i):
        tlast2 = m.Nfe.at(-2)
        return m.grad_u_L_ele[(tlast2, i)] == m.end_pen_grad_u_L[i]
    mod.grad_u_L_end = Constraint(mod.u_ind, rule=grad_u_L_end_rule)

    # Sum up the square of all penality vars
    def sum_square_end_pen_rule(m):
        term1 = sum(m.end_pen_x[i]**2 for i in m.tray)
        term2 = sum(m.end_pen_M[i]**2 for i in m.tray)
        term3 = sum(m.end_pen_grad_u_L[i]**2 for i in m.u_ind)

        return w_x*term1 + w_M*term2 + w_rg*term3
    # mod.sum_square_end_pen = Var(initialize=0.0, bounds=(0, None))
    # mod.con_sum_square_end_pen = Constraint(rule=sum_square_end_pen_rule)

    mod.sum_square_end_pen = Expression(rule=sum_square_end_pen_rule)


def construct_stabilizing_constraint(mod, stabilizing_delta, w_x, w_M, w_rg):
    mod.val_function = Var(initialize = 0.0)
    h = mod.h.value

    def con_val_function_rule(m):
        # w_x = 1.0
        # w_M = 1.0
        # w_rg = 1.0

        val_fun = sum(
            w_x*sum((m.x[(i+h, j)] - m.x[(i, j)])**2 for j in m.tray) +
            w_M*sum((m.M[(i+h, j)] - m.M[(i, j)])**2 for j in m.tray) +
            w_rg*sum((m.grad_u_L_ele[(i, k)])**2 for k in m.u_ind)
                for i in m.Nfe if i != m.Nfe.at(-1)
                        )

        return m.val_function == val_fun
    mod.con_val_function = Constraint(rule=con_val_function_rule)

    mod.pre_val_function = Param(initialize = 0.0, mutable = True)
    mod.pre_first_stage_cost = Param(initialize = 0.0, mutable = True)

    mod.stabilizing_delta = Param(initialize = stabilizing_delta)
    mod.stabilizing_con = Constraint(
        rule = lambda m:
            m.val_function - m.pre_val_function <= \
                -m.stabilizing_delta * m.pre_first_stage_cost
        )


def calculate_stage_cost_val_at_tp(m, tp, w_x, w_M, w_rg):
    h = m.h.value
    # w_x = 1.0
    # w_M = 1.0
    # w_rg = 1.0

    stage_cost_val = pyo_val(
        w_x*sum((m.x[(tp+h, j)] - m.x[(tp, j)])**2 for j in m.tray) +
        w_M*sum((m.M[(tp+h, j)] - m.M[(tp, j)])**2 for j in m.tray) +
        w_rg*sum((m.grad_u_L_ele[(tp, k)])**2 for k in m.u_ind)
                            )

    return stage_cost_val


def get_state_at_tp(mod, tp):
    states = dict()
    for sstr in state_strs:
        var = mod.find_component(sstr)
        for j in mod.tray:
            key = str(ComponentUID(var[:, j]))
            states[key] = var[tp, j].value

    D_cuid = ComponentUID(mod.D[:])
    states[str(D_cuid)] = mod.D[tp].value

    L1_cuid = ComponentUID(mod.L[:,1])
    states[str(L1_cuid)] = mod.L[tp,1].value

    return states


def get_control_at_tp(mod, tp):
    controls = dict()
    for cstr in control_strs:
        var = mod.find_component(cstr)
        key = str(ComponentUID(var[:]))
        controls[key] = var[tp].value

    return controls


def set_ics(mod, ic_state):
    for sstr in ["x", "M"]: # Do not use state_strs!!
        var = mod.find_component(sstr)
        for j in mod.tray:
            key = str(ComponentUID(var[:, j]))
            var[0, j].fix(ic_state[key])


def init_record_states(mod):
    record_state = dict()
    for sstr in state_strs:
        for j in mod.tray:
            var = mod.find_component(sstr)
            cuid = ComponentUID(var[:, j])
            record_state[str(cuid)] = []

    D_cuid = ComponentUID(mod.D[:])
    record_state[str(D_cuid)] = []

    L1_cuid = ComponentUID(mod.L[:,1])
    record_state[str(L1_cuid)] = []

    return record_state


# def save_states(record_state, mod, tp):
#     for sstr in state_strs:
#         for j in mod.tray:
#             var = mod.find_component(sstr)
#             key = str(ComponentUID(var[:, j]))
#             val = var[(tp, j)].value
#             record_state[key].append(val)

#     return record_state


def save_states(record_state, states):
    for key, val in record_state.items():
        val.append(states[key])

    return record_state


def init_record_ss_info(ss_info):
    # record_ss_info = copy.deepcopy(ss_info)
    record_ss_info = dict()
    for key, val in ss_info.items():
        record_ss_info[key] = []

    return record_ss_info


def save_ss_info(record_ss_info, ss_info):
    for key, val in record_ss_info.items():
        val.append(ss_info[key])

    return record_ss_info


def init_record_controls(mod):
    record_control = dict()
    for cstr in control_strs:
        var = mod.find_component(cstr)
        cuid = ComponentUID(var[:])
        record_control[str(cuid)] = []

    return record_control


def save_controls(record_control, controls):
    for key, val in record_control.items():
        val.append(controls[key])
    return record_control


# def init_record_econ_stage_cost(mod):
#     record_phi = []

#     return record_phi


def save_econ_stage_cost(record_phi, mod, tp):
    econ_val = pyo_val(mod.econ_stage_cost_expr[tp])
    print("economic stage cost: ", econ_val)
    record_phi.append(econ_val)

    return record_phi


def save_results(directory, res_states, res_controls, res_phi, res_ss_info):
    file_state = directory + "/res_states.json"
    with open(file_state, "w") as f1:
        json.dump(res_states, f1)

    file_control = directory + "/res_controls.json"
    with open(file_control, "w") as f2:
        json.dump(res_controls, f2)

    file_phi = directory + "/res_plant_phi.json"
    with open(file_phi, "w") as f3:
        json.dump({"phi":res_phi}, f3)

    file_ss = directory + "/res_ss_info.json"
    with open(file_ss, "w") as f4:
        json.dump(res_ss_info, f4)


def load_results(directory):
    file_state = directory + "/res_states.json"
    with open(file_state, "r") as f1:
        record_state = json.load(f1)

    file_control = directory + "/res_controls.json"
    with open(file_control, "r") as f2:
        record_control = json.load(f2)

    file_phi = directory + "/res_plant_phi.json"
    with open(file_phi, "r") as f3:
        res_phi = json.load(f3)
    record_plant_phi = res_phi["phi"]

    file_ss = directory + "/res_ss_info.json"
    with open(file_ss, "r") as f4:
        record_ss_info = json.load(f4)

    return record_state, record_control, record_plant_phi, record_ss_info


def solve_rg_variables_fixed_model_primals(mod, solver):
    for var in mod.component_objects(Var):
        try:
            var.fix()
        except AttributeError:
            var[...].fix()

    for con in mod.component_objects(Constraint):
        con.deactivate()

    mod.lamb_ele.unfix()
    mod.nu_ele.unfix()
    mod.eta_ele.unfix()
    mod.KKT1.activate()
    mod.KKT2.activate()
    mod.KKT4.activate()

    solve_problem("Initialize all lambda, nu, and eta", solver, mod)

    mod.grad_u_L_ele.unfix()
    mod.KKT3.activate()

    solve_problem("Initialize all grad_u_f", solver, mod)

    mod.val_function.unfix()
    mod.con_val_function.activate()

    solve_problem("Solve all rg expressions and value function", solver, mod)

    for var in mod.component_objects(Var):
        try:
            var.unfix()
        except AttributeError:
            var[...].unfix()

    for con in mod.component_objects(Constraint):
        con.activate()

    mod.x[0, :].fix()
    mod.M[0, :].fix()
    mod.Rec[0].fix()
    mod.Qr[0].fix()


def deactivate_rg_constraints(mod):
    mod.KKT1.deactivate()
    mod.KKT2.deactivate()
    mod.KKT3.deactivate()
    mod.KKT4.deactivate()

    mod.con_val_function.deactivate()
    mod.stabilizing_con.deactivate()


def save_value_function_info(dire,
                             record_value_function,
                             record_allowed_value_function,
                             record_x_sqsum, record_M_sqsum, record_rg_sqsum):
    file_val_function = dire + "/res_value_function.json"
    with open(file_val_function, "w") as f1:
        json.dump(record_value_function, f1)

    file_allowed_value_function = dire + "/res_allowed_value_function.json"
    with open(file_allowed_value_function, "w") as f2:
        json.dump(record_allowed_value_function, f2)

    file_x_sqsum = dire + "/res_x_sqsum.json"
    with open(file_x_sqsum, "w") as f3:
        json.dump(record_x_sqsum, f3)

    file_M_sqsum = dire + "/res_M_sqsum.json"
    with open(file_M_sqsum, "w") as f4:
        json.dump(record_M_sqsum, f4)

    file_rg_sqsum = dire + "/res_rg_sqsum.json"
    with open(file_rg_sqsum, "w") as f5:
        json.dump(record_rg_sqsum, f5)


def save_first_stage_cost_info(dire,
                               record_first_stage_cost,
                               record_x_sq_stagecost0, record_M_sq_stagecost0, record_rg_sq_stagecost0):
    file_first_stage_cost = dire + "/res_first_stage_cost.json"
    with open(file_first_stage_cost, "w") as f1:
        json.dump(record_first_stage_cost, f1)

    file_x_sq_stagecost0 = dire + "/res_x_sq_stagecost0.json"
    with open(file_x_sq_stagecost0, "w") as f2:
        json.dump(record_x_sq_stagecost0, f2)

    file_M_sq_stagecost0 = dire + "/res_M_sq_stagecost0.json"
    with open(file_M_sq_stagecost0, "w") as f3:
        json.dump(record_M_sq_stagecost0, f3)

    file_rg_sq_stagecost0 = dire + "/res_rg_sq_stagecost0.json"
    with open(file_rg_sq_stagecost0, "w") as f4:
        json.dump(record_rg_sq_stagecost0, f4)
