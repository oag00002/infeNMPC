import pyomo.environ as pyo

def econ_stage_cost1(mod, fe_set, dynamic=False):
    # original eco5
    c1 = 0.5
    p1 = 1.5
    p2 = 1.0

    def econ_stage_cost_rule(m, i):
        # if dynamic:
        #     if i == fe_set.at(-1):
        #         return pyo.Constraint.Skip
        if dynamic:
            h = mod.h.value
            if i == mod.Nfe.last():
                return pyo.Constraint.Skip
            else:
                return (c1*mod.Qr[i+h]) - (p1*mod.D[i] + p2*mod.L[i,1])
        else:
            return (c1*mod.Qr[i]) - (p1*mod.D[i] + p2*mod.L[i,1])
    return pyo.Expression(fe_set, rule=econ_stage_cost_rule)

def econ_bound1(mod, dynamic=False, relaxBounds=False):
    x42lb = 0.98
    x1ub = 0.02

    if dynamic:
        if relaxBounds:
            x42_lbVar = pyo.Var(mod.t, initialize=0.0, bounds=(0, None))
            def x42_relax_rule(m,i):
                if i==0:
                    return pyo.Constraint.Skip
                else:
                    return m.x[i,42] + m.x42_lbVar[i] >= x42lb
            x42_lbCon = pyo.Constraint(mod.t, rule=x42_relax_rule)

            x1_ubVar = pyo.Var(mod.t, initialize=0.0, bounds=(0, None))
            def x1_relax_rule(m,i):
                if i==0:
                    return pyo.Constraint.Skip
                else:
                    return m.x[i,1] <= 0.02 + m.x1_ubVar[i]
            x1_ubCon = pyo.Constraint(mod.t, rule=x1_relax_rule)

            return x42_lbVar, x42_lbCon, x1_ubVar, x1_ubCon

        else:
            mod.x[:, 42].setlb(x42lb)
            mod.x[:, 1].setub(x1ub)
            mod.x[0, 42].setlb(0.0)
            mod.x[0, 1].setub(1.0)
    else:
        mod.x[:, 42].setlb(x42lb)
        mod.x[:, 1].setub(x1ub)

def tighter_control_bound1_for_eNMPC(mod):
    RecUb = 3 #10.0#5.0
    QrUb = 20 #25.0

    Rec_ubVar = pyo.Var(mod.t, initialize=0.0, bounds=(0, None))
    def Rec_relax_rule(m,i):
        if i==0:
            return pyo.Constraint.Skip
        else:
            return m.Rec[i] <= RecUb + m.Rec_ubVar[i]
    Rec_ubCon = pyo.Constraint(mod.t, rule=Rec_relax_rule)

    Qr_ubVar = pyo.Var(mod.t, initialize=0.0, bounds=(0, None))
    def Qr_relax_rule(m,i):
        if i==0:
            return pyo.Constraint.Skip
        else:
            return m.Qr[i] <= QrUb + m.Qr_ubVar[i]
    Qr_ubCon = pyo.Constraint(mod.t, rule=Qr_relax_rule)

    return Rec_ubVar, Rec_ubCon, Qr_ubVar, Qr_ubCon


def calibrated_bound_expression_for_KKT_conditions(mod):
    mu_x42 = 1.0E-3
    mu_x1 = 1.0E-3

    def h_x42_rule(m, i):
        if i == mod.Nfe.at(-1):
            return pyo.Constraint.Skip
        else:
            # return 0.98 -m.x[i, 42] - 9.58956127151866e-11 # for mu=1.0E-8.6
            return 0.98 -mod.x[i, 42] - 3.8176731069480156e-05 # for mu=1.0E-3
    h_x42 = pyo.Expression(mod.Nfe, rule=h_x42_rule)

    def h_x1_rule(m, i):
        if i == mod.Nfe.at(-1):
            return pyo.Constraint.Skip
        else:
            # return m.x[i, 1] - 0.02 - 7.931090414593066e-09 # for mu=1.0E-8.6
            return mod.x[i, 1] - 0.02 - 0.0031574239662551414 # for mu=1.0E-3
    h_x1 = pyo.Expression(mod.Nfe, rule=h_x1_rule)

    return mu_x42, h_x42, mu_x1, h_x1

#-----------------------------------------------------------------------------

cost_dict = {1:econ_stage_cost1}
def dist_econ_stage_cost(cost_num):
    return cost_dict[cost_num]

bound_dict = {1:econ_bound1}
def dist_econ_bound(cost_num):
    return bound_dict[cost_num]

tighter_control_bound_dict = {1:tighter_control_bound1_for_eNMPC}
def dist_tighter_control_bound_for_eNMPC(cost_num):
    return tighter_control_bound_dict[cost_num]