"""
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
        return m.dcadt[i] == m.q[i]/m.V*(m.caf[i] - m.ca[i]) - m.k*m.ca[i]

def _cbdot(m,i):
    if i == 0:
        return Constraint.Skip
    else:
        return m.dcbdt[i] == m.q[i]/m.V*(-m.cb[i]) + m.k*m.ca[i]

def CSTR_daemodel(time_period):

    m = ConcreteModel()
    m.h = Param(initialize = 1.0)
    m.t = ContinuousSet(bounds = time_period)

    #Parameters
    m.V = Param(initialize = 10.0)
    m.k = Param(initialize = 1.2)
    m.caf = Param(m.t, initialize = 1.0, mutable = True)

    #States
    m.ca = Var(m.t, initialize = 0.5, bounds = (0, None))
    m.cb = Var(m.t, initialize = 0.5, bounds = (0, None))
    m.dcadt = DerivativeVar(m.ca, wrt = m.t)
    m.dcbdt = DerivativeVar(m.cb, wrt = m.t)

    #Controls
    m.q = Var(m.t, initialize = 12., bounds = (0, None))

    #Constraints
    m.cadot = Constraint(m.t, rule = _cadot)
    m.cbdot = Constraint(m.t, rule = _cbdot)

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

def rg_cadot(m,i):
    h = m.Nfe.at(2) - m.Nfe.at(1)
    if i == m.Nfe.at(-1):
        return Constraint.Skip
    else:
        return m.ca[i] - (m.q[i+h]/m.V*(m.caf[i] - m.ca[i]) - m.k*m.ca[i])

def rg_cbdot(m,i):
    h = m.Nfe.at(2) - m.Nfe.at(1)
    if i == m.Nfe.at(-1):
        return Constraint.Skip
    else:
        return m.cb[i] - (m.q[i+h]/m.V*(-m.cb[i]) + m.k*m.ca[i])

def set_caf_after_discretize(mod, caf_val):
    for ind in mod.caf:
        mod.caf[ind] = caf_val

