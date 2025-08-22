from workspace.models.diehl_dist_model import dist_daemodel, create_model_hard_bounds
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

mod = dist_daemodel(nfe = 1)#, delete_unused_var = True)
augment_de_to_ss_eq(mod)
create_model_hard_bounds(mod)
# Tracking stage cost
# mod.obj = Objective(expr = (0.48289976241629246*(mod.T[60,29] - 345.22)**2 +
#                             0.2032811003140813*(mod.T[60,14] - 356.23)**2))
# mod.obj1 = Objective(expr = (0.48289976241629246*(mod.T[60,29] - 346.5)**2 +
#                               0.2032811003140813*(mod.T[60,14] - 354.6)**2))

# mod.obj_e = Objective(expr =
#                       (mod.D[60]-17.7)**2 +
#                         100000*(mod.x[60, 42]-0.980)**2 +
#                       (mod.Qr[60]-19.0)**2)

# mod.obj_e2 = Objective(expr =
#                       (mod.D[60]-18.76)**2)
# mod.x[60, 42].fix(0.98)

# mod.x[60, 42].setlb(0.98)

# economic stage cost
# mod.obj2_scaled = Objective(expr = -5.0E1*mod.x[60, 42] - 1.0E0*mod.D[60] + mod.Qr[60])
# obj2 fails, I should add penality on the purity requirement

# mod.obj3 = Objective(expr= -1.0E-1*mod.D[60] + 1.0E-1*mod.Qr[60])
# mod.x[60, 42].setlb(0.98)
# construct_stage_cost(mod, mod.t, dynamic=False)
# mod.obj4 = Objective(expr=mod.econ_stage_cost_expr[60])

# mod.obj5 = Objective(expr = -mod.D[60])
# mod.Rec[60].setub(3.0)
# mod.D[60].setlb(18.0)

# mod.obj6 = Objective(expr = -mod.D[60])
# mod.Rec[60].setub(1.5)
# mod.x[60, 42].setlb(0.99)

# mod.obj_e2 = Objective(expr =
#                       (mod.D[60]-18.54)**2)
# mod.x[60, 42].fix(0.99)

# # eco4
# # mod.obj = Objective(expr=mod.Qr[60])
# mod.obj = Objective(expr = (
#     (mod.D[60]-18.0)**2 + 1.0E5*(mod.x[60,42]-0.98)**2 + (mod.Qr[60]-17.25)**2
#     ))

# eco5
c1 = 0.5
p1 = 1.5
p2 = 1.0

mod.x[60,42].setlb(0.98)
mod.x[60,1].setub(0.02)

# mod.obj_eco5 = Objective(
#     expr=(c1*mod.Qr[60]) - (p1*mod.D[60] + p2*mod.L[60,1])
#     )
mod.obj_ic = Objective(
    expr=1.0E2*(mod.x[60,42]-0.982)**2 + 1.0E2*(mod.x[60,1]-0.198)**2 + (mod.Qr[60]-17.16)**2
    )

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

for obj in mod.component_objects(Objective):
    if obj.active:
        print("economic_cost: ", pyo_val(obj))


# # Export json file
# scalar_var, dae_var = flatten_dae_components(mod,
#                                              mod.t,
#                                              ctype = Var)

# data_dict = dict()
# for var in dae_var:
#     cuid = ComponentUID(var.referent)
#     val = var[60].value
#     data_dict[str(cuid)] = val

# import json
# with open("dist_ss1_scaled.json", "w") as f1:
#     json.dump(data_dict, f1)
