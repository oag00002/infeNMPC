from workspace.models.diehl_dist_model_w_t0_con_feed_as_u import (
    dist_daemodel, create_model_hard_bounds
)
from workspace.dynamics import augment_de_to_ss_eq
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    )
from pyomo.environ import (
    Constraint,
    Expression,
    Objective,
    SolverFactory,
    Var,
    ComponentUID,
    Reference,
    log,
    Param,
    Suffix,
    )
from pyomo.environ import value as pyo_val
from pyomo.dae import (
    DerivativeVar
    )
from pyomo.dae.flatten import flatten_dae_components

mod = dist_daemodel(nfe=1)#, delete_unused_var = True)
augment_de_to_ss_eq(mod)
create_model_hard_bounds(mod)
t0 = mod.t.first()
scalar_con, dae_con = flatten_dae_components(mod,
                                             mod.t,
                                             ctype = Constraint)
for con in dae_con:
    con[t0].deactivate()

# Tracking stage cost
# mod.obj = Objective(expr = (0.48289976241629246*(mod.T[60,29] - 345.22)**2 +
#                             0.2032811003140813*(mod.T[60,14] - 356.23)**2))
# mod.obj1 = Objective(expr = (0.48289976241629246*(mod.T[60,29] - 346.5)**2 +
#                               0.2032811003140813*(mod.T[60,14] - 354.6)**2))

# # eco3
# mod.obj_eco3 = Objective(expr = -mod.D[60])
# mod.feed.fix()
# mod.x[60, 42].setlb(0.98)
# mod.Rec[60].setub(1.5)

# mod.obj_e2 = Objective(expr =
#                       (mod.D[60]-18.54)**2)
# mod.x[60, 42].fix(0.99)

# # eco4
# mod.obj_eco4 = Objective(expr=mod.Qr[60])
# mod.feed.fix()
# mod.x[60, 42].setlb(0.98)
# mod.D[60].setlb(18.0)

# eco5 (not use)
c1 = 1.0
c2 = 0.5
p1 = 1.5
p2 = 1.0
mod.obj_eco5 = Objective(expr=
                          # (c1*mod.feed[60,21] + c2*mod.Qr[60]) -
                          (c2*mod.Qr[60]) -
                          (p1*mod.D[60] + p2*mod.L[60,1])
                          )
# mod.Rec.fix(1.24)
# mod.feed.setub(100.0)
mod.feed[60,21].fix(57.5294)
mod.x[60, 42].setlb(0.98)
mod.x[60, 1].setub(0.02)

# # eco6 (not use)
# mod.obj_eco6 = Objective(expr=mod.Qr[60])
# mod.x[60, 42].setlb(0.98)
# mod.DF_con = Constraint(expr= mod.D[60] >= 0.3*mod.feed[60,21])


### Declare all suffixes
# Ipopt bound multipliers (obtained from solution)
mod.ipopt_zL_out = Suffix(direction=Suffix.IMPORT)
mod.ipopt_zU_out = Suffix(direction=Suffix.IMPORT)
# Ipopt bound multipliers (sent to solver)
mod.ipopt_zL_in = Suffix(direction=Suffix.EXPORT)
mod.ipopt_zU_in = Suffix(direction=Suffix.EXPORT)
# Obtain dual solutions from first solve and send to warm start
mod.dual = Suffix(direction=Suffix.IMPORT_EXPORT)


dof = degrees_of_freedom(mod)
print("dof = ", dof)

solver = SolverFactory("ipopt")
# solver.options["halt_on_ampl_error"] = "yes"
# solver.options["max_iter"] = 14
solver.solve(mod, tee = True)

print(mod.T[60, 29].value)
print(mod.T[60, 14].value)

mod.Rec[60].pprint()
mod.Qr[60].pprint()
mod.x[60, 42].pprint()
mod.x[60, 1].pprint()
mod.D[60].pprint()
mod.L[60,1].pprint()
mod.feed[60,21].pprint()

for obj in mod.component_objects(Objective):
    if obj.active:
        print("economic_cost: ", pyo_val(obj))

# # Export json file
# scalar_var, dae_var = flatten_dae_components(mod,
#                                               mod.t,
#                                               ctype = Var)
# data_dict = dict()
# for var in dae_var:
#     cuid = ComponentUID(var.referent)
#     val = var[60].value
#     data_dict[str(cuid)] = val

# import json
# with open("dist_ss1_scaled.json", "w") as f1:
#     json.dump(data_dict, f1)
