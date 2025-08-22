"""
Note that this is a DAE model!!!

CSTR with reaction:
    A ----> B, r = k*ca
from "A Lyapunov Function for Economic Optimizing Model Predictive Control" by Diehl et al.

@author: Kuan-Han Lin
"""

from pyomo.environ import (
    ConcreteModel,
    Param,
    Var,
    Constraint,
    Expression,
    )
from pyomo.dae import (
    ContinuousSet,
    DerivativeVar,
    )

#Main model--------------------------------------------------------------------
#Equations
def _cadot(m,i):
    if i == 0:
        return Constraint.Skip
    else:
        return m.dcadt[i] == m.caa[i] - m.k*m.ca[i]

def _cbdot(m,i):
    if i == 0:
        return Constraint.Skip
    else:
        return m.dcbdt[i] == m.cbb[i] + m.k*m.ca[i]

def _caa_con(m,i):
    if i == 0:
        return Constraint.Skip
    else:
        return m.q[i]/m.V*(m.caf[i] - m.ca[i]) - m.caa[i] == 0.0

def _cbb_con(m,i):
    if i == 0:
        return Constraint.Skip
    else:
        return m.q[i]/m.V*(-m.cb[i]) - m.cbb[i] == 0.0


def CSTR_daemodel(time_period):

    m = ConcreteModel()
    m.h = Param(initialize = 1.0)
    m.t = ContinuousSet(bounds = time_period)

    #Parameters
    m.V = Param(initialize = 10.0)
    m.k = Param(initialize = 1.2)
    # m.caf = Param(m.t, initialize = 1.0, mutable = True)
    m.caf = Var(m.t, initialize = 1.0)
    m.caf.fix()


    #States
    m.ca = Var(m.t, initialize = 0.5, bounds = (0, None))
    m.cb = Var(m.t, initialize = 0.5, bounds = (0, None))
    m.dcadt = DerivativeVar(m.ca, wrt = m.t)
    m.dcbdt = DerivativeVar(m.cb, wrt = m.t)

    #Algebraic variables
    m.caa = Var(m.t, initialize = 0.6)
    m.cbb = Var(m.t, initialize = -0.6)

    #Controls
    m.q = Var(m.t, initialize = 12., bounds = (0, None))

    #Constraints
    m.cadot = Constraint(m.t, rule = _cadot)
    m.cbdot = Constraint(m.t, rule = _cbdot)
    m.caa_con = Constraint(m.t, rule = _caa_con)
    m.cbb_con = Constraint(m.t, rule = _cbb_con)

    states = [m.ca, m.cb]
    controls = [m.q]

    return m, states, controls

#create bound constraints------------------------------------------------------
def create_soft_bound_constraint(model, var, bound, val):
    rvar_name = "r" + var.name + bound
    model.add_component(rvar_name, Var(model.Nfe,
                                       initialize = 0.0,
                                       bounds = (0., None)))
    rvar = model.find_component(rvar_name)

    if bound == "lb":
        # def _rvar_lb(model, i):
        #     if i == 0:
        #         return Constraint.Skip
        #     else:
        #         return bound_value <= item[i] + rvar[i]
        rcon_name = "soft_" + var.name + "_" + bound
        model.add_component(rcon_name, Constraint(
            model.Nfe,
            rule = lambda model, i: Constraint.Skip if i==0
                                    else val <= var[i] + rvar[i])
                            )
        # del _rvar_lb

    elif bound == "ub":
        # def _rvar_ub(model, i):
        #     if i == 0:
        #         return Constraint.Skip
        #     else:
        #         return bound_value + rvar[i] >= item[i]
        rcon_name = "soft_" + var.name + "_" + bound
        model.add_component(rcon_name, Constraint(
            model.Nfe,
            rule = lambda model, i: Constraint.Skip if i==0
                                    else val + rvar[i] >= var[i])
                            )
        # del _rvar_ub
    else:
        raise RuntimeError("Bound can only be either lb or ub!")

    rcon = model.find_component(rcon_name)

    return rvar, rcon

def create_model_bounds(model, bound_type):
    cb_bounds = {"lb": 0.45, "ub": 1.0}
    q_bounds = {"lb": 10.0, "ub": 20.0}
    var_w_bounds = {model.cb: cb_bounds,
                    model.q: q_bounds}

    if bound_type == "hard":
        for var, bounds in var_w_bounds.items():
            var.setlb(bounds["lb"])
            var.setub(bounds["ub"])

    elif bound_type == "soft":
        rvar_set = set()
        rcon_set = set()
        for var, bounds in var_w_bounds.items():
            for bound, val in bounds.items():
                rvar, rcon = create_soft_bound_constraint(model,
                                                          var,
                                                          bound,
                                                          val)
                rvar_set.add(rvar)
                rcon_set.add(rcon)

        return rvar_set, rcon_set

    else:
        raise RuntimeError("unrecognizable bound type!")

def get_sample_time(blo):
    h = blo.find_component("h")
    if h is None:
        h = blo.parent_block().find_component("h")
    return h.value

def get_nfe_set(blo):
    Nfe = blo.find_component("Nfe")
    if Nfe is None:
        Nfe = blo.parent_block().find_component("Nfe")
    return Nfe

def rg_cadot(b,i):
    h = get_sample_time(b)
    Nfe = get_nfe_set(b)
    m = b.parent_block()
    if i == Nfe.at(-1):
        return Constraint.Skip
    else:
        return m.caa[i] - m.k*m.ca[i]

def rg_cbdot(b,i):
    h = get_sample_time(b)
    Nfe = get_nfe_set(b)
    m = b.parent_block()
    if i == Nfe.at(-1):
        return Constraint.Skip
    else:
        return m.cbb[i] + m.k*m.ca[i]

def rg_caa_con(b,i):
    h = get_sample_time(b)
    Nfe = get_nfe_set(b)
    m = b.parent_block()
    if i == Nfe.at(-1):
        return Constraint.Skip
    else:
        return m.q[i+h]/m.V*(m.caf[i] - m.ca[i]) - m.caa[i]

def rg_cbb_con(b,i):
    h = get_sample_time(b)
    Nfe = get_nfe_set(b)
    m = b.parent_block()
    if i == Nfe.at(-1):
        return Constraint.Skip
    else:
        return m.q[i+h]/m.V*(-m.cb[i]) - m.cbb[i]

def construct_rg_model_expr(blo):
    par = blo.parent_block()
    blo.rg_cadot = Expression(par.Nfe, rule = rg_cadot)
    blo.rg_cbdot = Expression(par.Nfe, rule = rg_cbdot)
    blo.rg_caa_con = Expression(par.Nfe, rule = rg_caa_con)
    blo.rg_cbb_con = Expression(par.Nfe, rule = rg_cbb_con)

def set_caf_after_discretize(mod, caf_val):
    for ind in mod.caf:
        mod.caf[ind] = caf_val
