from pyomo.environ import (
    Set,
    Var,
    Constraint,
    TransformationFactory,
    Param,
    Objective,
    Expression,
    ComponentUID,
    )
from pyomo.environ import value as pyo_val
from pyomo.dae import (
    DerivativeVar,
    ContinuousSet,
    )
from pyomo.dae.flatten import flatten_dae_components

def augment_de_to_ss_eq(model):
    for dv in model.component_objects(DerivativeVar):
        for index, item in dv.items():
            item.fix(0.0)


def discretize_and_create_fe_set(model, nfe, ncp,
                                 is_control_model,
                                 control_vars = None):
    discretizer = TransformationFactory("dae.collocation")
    discretizer.apply_to(model,
                         nfe = nfe,
                         ncp = ncp,
                         scheme = "LAGRANGE-RADAU")
    if is_control_model:
        for cvar in control_vars:
            discretizer.reduce_collocation_points(model,
                                                  var = cvar,
                                                  ncp = 1,
                                                  contset = model.t)
    steps = model.t.get_finite_elements()
    model.Nfe = ContinuousSet(initialize = steps)
    # model.h = Param(initialize = model.Nfe.at(2) - model.Nfe.at(1))


def shift_dae_vars_by_dt(dae_vars, t_set, dt):
    seen = set()
    t0 = t_set.first()
    tf = t_set.last()
    for var in dae_vars:
        if id(var[t0]) in seen:
            continue
        else:
            seen.add(id(var[t0]))
        for t in t_set:
            ts = t + dt
            idx = t_set.find_nearest_index(ts)
            if idx is None:
                # ts is outside the controller's horizon
                var[t].set_value(var[tf].value)
            else:
                ts = t_set.at(idx)
                var[t].set_value(var[ts].value)

def construct_tracking_stage_cost(
        t_set, weights, variables, setpoints):
    if len(weights) != len(variables):
        raise RuntimeError("num of weights and num of variables do not match!")
    if len(setpoints) != len(variables):
        raise RuntimeError("num of setpoints and num of variables do not match!")

    def tracking_stage_cost(m, i):
        return sum(weight*(variable[i]-setpoint)**2
                   for weight, variable, setpoint in zip(weights, variables, setpoints))

    return Expression(t_set, rule=tracking_stage_cost)

def get_cuid_of_var(var):
    if var.is_reference():
        return ComponentUID(var.referent)
    else:
        return ComponentUID(var)

def set_initial_conds(t0, states, ic_info):
    for state in states:
        cuid = get_cuid_of_var(state)
        state[t0].fix(ic_info[str(cuid)])

def initialize_vars(mod, var_info):
    for cuid_str, val in var_info.items():
        cuid = ComponentUID(cuid_str)
        var = cuid.find_component_on(mod)
        var[:].set_value(val)

def get_comp_value_at_tp(components, tp):
    comp_dict = dict()
    for comp in components:
        cuid = get_cuid_of_var(comp)
        comp_dict[str(cuid)] = pyo_val(comp[tp])
    return comp_dict

def construct_L2relaxed_endpoint_cons(states, steadystates, tf):
    if len(states) != len(steadystates):
        raise RuntimeError("Num of states and steady-states do not match!")
    endpt_set = Set(initialize=[i for i in range(len(states))])
    endpt_var = Var(endpt_set, initialize=0.0)
    def endpoint_con_rule(m, i):
        return states[i][tf] - steadystates[i] == endpt_var[i]
    endpt_con = Constraint(endpt_set, rule=endpoint_con_rule)

    return endpt_set, endpt_var, endpt_con

def construct_L1relaxed_endpoint_cons(states, steadystates, tf):
    if len(states) != len(steadystates):
        raise RuntimeError("Num of states and steady-states do not match!")
    endpt_set = Set(initialize=[i for i in range(len(states))])
    endpt_varType = Set(initialize=["p", "n"])
    endpt_var = Var(endpt_set, endpt_varType, initialize=0.0, bounds=(0.0, None))
    def endpoint_con_rule(m, i):
        return states[i][tf] - steadystates[i] == endpt_var[i,"p"] - endpt_var[i,"n"]
    endpt_con = Constraint(endpt_set, rule=endpoint_con_rule)

    return endpt_set, endpt_varType, endpt_var, endpt_con

def retrieve_dae_vars_at_tp(mod, tp):
    scalar_vars, dae_vars = flatten_dae_components(mod, mod.t, Var)
    var_dict = dict()
    for var in dae_vars:
        cuid = ComponentUID(var.referent)
        var_dict[str(cuid)] = var[tp].value
    return var_dict