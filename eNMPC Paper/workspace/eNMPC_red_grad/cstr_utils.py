from pyomo.environ import (
    Var,
    Objective,
    SolverFactory,
    )
from pyomo.core.base.var import IndexedVar
from workspace.dynamics import augment_de_to_ss_eq
from workspace.models.diehl_cstr_model import (
    CSTR_daemodel,
    create_model_bounds,
    )


def obtain_optimal_ss(caf = 1.):
    time_period = (0, 1)
    e, states, controls = CSTR_daemodel(time_period)
    augment_de_to_ss_eq(e)

    create_model_bounds(e, bound_type = "hard")
    for ind, item in e.caf.items():
        item.set_value(caf)
    e.obj_ss = Objective(expr = -e.q[1]*(2*e.cb[1]-0.5))

    solver = SolverFactory("ipopt")
    print("Solve steady-state problem")
    results = solver.solve(e, tee = True)
    print("optiaml steady state")
    print("ca: ", e.ca[1].value, " cb: ", e.cb[1].value)
    print("q: ", e.q[1].value)

    ss = {"ca": e.ca[1].value,
          "cb": e.cb[1].value,
          "q": e.q[1].value,
          "Lss": e.obj_ss.expr()}

    return ss


def set_caf(model, caf_val):
    for index, item in model.caf.items():
        item.set_value(caf_val)
    if type(model.caf) is IndexedVar:
        model.caf.fix()

def save_state(record, data):
    for i in record.keys():
        record[i].append(data[i])
    return record


def save_ss_state(record, ss_data):
    for i in record.keys():
        record[i].append(ss_data[i])
    return record


def set_ics(model, dict_ics):
    model.ca[0].fix(dict_ics["ca"])
    model.cb[0].fix(dict_ics["cb"])


