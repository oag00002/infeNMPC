from pyomo.environ import (
    Param,
    Var,
    Constraint,
    ComponentUID,
    )
from pyomo.environ import value as pyo_val
import json

def construct_eNMPC_sc_stabilizing_constrint(
        mod, ss, delta, w_x, w_M, w_Rec, w_Qr
    ):
    mod.val_function = Var(initialize = 0.0)
    h = mod.h.value

    def con_val_function_rule(m):
        val_fun = \
        sum(
            w_x*sum((m.x[(i, j)] - ss[str(ComponentUID(mod.x[:,j]))])**2 for j in m.tray) +
            w_M*sum((m.M[(i, j)] - ss[str(ComponentUID(mod.M[:,j]))])**2 for j in m.tray)
                for i in m.Nfe #if i != m.Nfe.at(-1)
                        ) + \
        sum(
            w_Rec*(m.Rec[i+h] - ss[str(ComponentUID(m.Rec[:]))])**2 +
            w_Qr*(m.Qr[i+h] - ss[str(ComponentUID(m.Qr[:]))])**2
                for i in m.Nfe if i != m.Nfe.at(-1)
                        )
        return m.val_function == val_fun
    mod.con_val_function = Constraint(rule=con_val_function_rule)

    mod.pre_val_function = Param(initialize = 0.0, mutable = True)
    mod.pre_first_stage_cost = Param(initialize = 0.0, mutable = True)

    mod.stabilizing_delta = Param(initialize = delta)
    mod.stabilizing_con = Constraint(
        rule = lambda m:
            m.val_function - m.pre_val_function <= \
                -m.stabilizing_delta * m.pre_first_stage_cost
        )

def construct_eNMPC_sc_relaxed_stabilizing_constrint(
        mod, ss, delta, w_x, w_M, w_Rec, w_Qr
    ):
    mod.val_function = Var(initialize = 0.0)
    h = mod.h.value

    def con_val_function_rule(m):
        val_fun = \
        sum(
            w_x*sum((m.x[(i, j)] - ss[str(ComponentUID(mod.x[:,j]))])**2 for j in m.tray) +
            w_M*sum((m.M[(i, j)] - ss[str(ComponentUID(mod.M[:,j]))])**2 for j in m.tray)
                for i in m.Nfe #if i != m.Nfe.at(-1)
                        ) + \
        sum(
            w_Rec*(m.Rec[i+h] - ss[str(ComponentUID(m.Rec[:]))])**2 +
            w_Qr*(m.Qr[i+h] - ss[str(ComponentUID(m.Qr[:]))])**2
                for i in m.Nfe if i != m.Nfe.at(-1)
                        )
        return m.val_function == val_fun
    mod.con_val_function = Constraint(rule=con_val_function_rule)

    mod.pre_val_function = Param(initialize = 0.0, mutable = True)
    mod.pre_first_stage_cost = Param(initialize = 0.0, mutable = True)
    mod.relaxed_sc_var = Var(initialize = 0.0, bounds=(0.0, None))

    mod.stabilizing_delta = Param(initialize = delta)
    mod.stabilizing_con = Constraint(
        rule = lambda m:
            m.val_function - m.pre_val_function <= \
                -m.stabilizing_delta * m.pre_first_stage_cost + m.relaxed_sc_var
        )

def calculate_eNMPC_sc_stage_cost_val_at_tp(
        m, ss, tp, w_x, w_M, w_Rec, w_Qr
    ):
    h = m.h.value

    stage_cost_val = pyo_val(
        w_x*sum((m.x[(tp, j)] - ss[str(ComponentUID(m.x[:,j]))])**2 for j in m.tray) +
        w_M*sum((m.M[(tp, j)] - ss[str(ComponentUID(m.M[:,j]))])**2 for j in m.tray) +
        w_Rec*(m.Rec[tp+h] - ss[str(ComponentUID(m.Rec[:]))])**2 +
        w_Qr*(m.Qr[tp+h] - ss[str(ComponentUID(m.Qr[:]))])**2
                            )

    return stage_cost_val

def save_value_function_info(dire,
                             record_value_function,
                             record_allowed_value_function,
                             record_x_sqsum, record_M_sqsum,
                             record_Rec_sqsum, record_Qr_sqsum):
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

    file_Rec_sqsum = dire + "/res_Rec_sqsum.json"
    with open(file_Rec_sqsum, "w") as f5:
        json.dump(record_Rec_sqsum, f5)

    file_Qr_sqsum = dire + "/res_Qr_sqsum.json"
    with open(file_Qr_sqsum, "w") as f6:
        json.dump(record_Qr_sqsum, f6)


def save_first_stage_cost_info(dire,
                               record_first_stage_cost,
                               record_x_sq_stagecost0, record_M_sq_stagecost0,
                               record_Rec_sq_stagecost0, record_Qr_sq_stagecost0):
    file_first_stage_cost = dire + "/res_first_stage_cost.json"
    with open(file_first_stage_cost, "w") as f1:
        json.dump(record_first_stage_cost, f1)

    file_x_sq_stagecost0 = dire + "/res_x_sq_stagecost0.json"
    with open(file_x_sq_stagecost0, "w") as f2:
        json.dump(record_x_sq_stagecost0, f2)

    file_M_sq_stagecost0 = dire + "/res_M_sq_stagecost0.json"
    with open(file_M_sq_stagecost0, "w") as f3:
        json.dump(record_M_sq_stagecost0, f3)

    file_Rec_sq_stagecost0 = dire + "/res_Rec_sq_stagecost0.json"
    with open(file_Rec_sq_stagecost0, "w") as f4:
        json.dump(record_Rec_sq_stagecost0, f4)

    file_Qr_sq_stagecost0 = dire + "/res_Qr_sq_stagecost0.json"
    with open(file_Qr_sq_stagecost0, "w") as f4:
        json.dump(record_Qr_sq_stagecost0, f4)