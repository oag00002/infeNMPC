import pyomo.environ as pyo

# # eco5
# c1 = 0.5
# p1 = 1.5
# p2 = 1.0
# def dist_econ_stage_cost(mod, fe_set, dynamic,):
#     h = mod.h.value

#     def econ_stage_cost_rule(m, i):
#         # if dynamic:
#         #     if i == fe_set.at(-1):
#         #         return pyo.Constraint.Skip
#         if dynamic:
#             if i == mod.Nfe.last():
#                 return pyo.Constraint.Skip
#             else:
#                 return (c1*mod.Qr[i+h]) - (p1*mod.D[i] + p2*mod.L[i,1])
#         else:
#             return (c1*mod.Qr[i]) - (p1*mod.D[i] + p2*mod.L[i,1])
#     mod.econ_stage_cost_expr = pyo.Expression(fe_set,
#                                               rule=econ_stage_cost_rule)

#------------------------------------------------------------------------------
# # eco4
# def dist_econ_stage_cost(mod, fe_set, dynamic,):
#     h = mod.h.value

#     def econ_stage_cost_rule(m, i):
#         # if dynamic:
#         #     if i == fe_set.at(-1):
#         #         return pyo.Constraint.Skip
#         if dynamic:
#             if i == mod.Nfe.last():
#                 return pyo.Constraint.Skip
#             else:
#                 return mod.Qr[i+h]
#         else:
#             return mod.Qr[i]
#     mod.econ_stage_cost_expr = pyo.Expression(fe_set,
#                                               rule=econ_stage_cost_rule)


# -----------------------------------------------------------------------------
#eco3
def state_cons_bc_stage_cost(mod):
    mod.x[:, 42].setlb(0.98)
    # mod.x[:, 42].setlb(0.985)
    mod.Rec[:].setub(1.5)
    # print("No specific bound constraints now!")


def dist_econ_stage_cost(mod, fe_set, dynamic,):
    h = mod.h.value

    def econ_stage_cost_rule(m, i):
        # if dynamic:
        #     if i == fe_set.at(-1):
        #         return pyo.Constraint.Skip

        return -mod.D[i]

    mod.econ_stage_cost_expr = pyo.Expression(fe_set,
                                              rule=econ_stage_cost_rule)


def construct_stage_cost(mod, fe_set, dynamic):
    # state_cons_bc_stage_cost(mod)
    dist_econ_stage_cost(mod, fe_set, dynamic)


def construct_rg_stage_cost(mod, fe_set):#, mu):
    h = mod.h.value
    # mod.mu = pyo.Param(initialize = mu)

    def rg_stage_cost_rule(m, i):
        if i == fe_set.at(-1):
            return pyo.Constraint.Skip
        else:
            return -mod.D[i]
    mod.rg_econ_stage_cost = pyo.Expression(fe_set, rule=rg_stage_cost_rule)

# -----------------------------------------------------------------------------
# D_coff = 1.0E-1
# Qr_coff = 1.0E-1

# # def state_cons_bc_stage_cost(mod):
# #     mod.x[:, 42].setlb(0.98)
# #     # mod.D[:].setlb(17.0+1.0E-10)
# #     # print("No specific bound constraints now!")


# def dist_econ_stage_cost(mod, fe_set, dynamic,):
#     h = mod.h.value

#     def econ_stage_cost_rule(m, i):
#         if dynamic:
#             if i == fe_set.at(-1):
#                 return pyo.Constraint.Skip
#             input_index = i+h
#         else:
#             input_index = i

#         return (-D_coff*mod.D[i] +
#                 Qr_coff*mod.Qr[input_index])# -
#                 # mod.mu*pyo.log(mod.x[i, 42]-0.98))

#     mod.econ_stage_cost_expr = pyo.Expression(fe_set,
#                                               rule=econ_stage_cost_rule)


# # def construct_stage_cost(mod, fe_set, dynamic):
# #     state_cons_bc_stage_cost(mod)
# #     dist_econ_stage_cost(mod, fe_set, dynamic)


# def construct_rg_stage_cost(mod, fe_set, mu):
#     h = mod.h.value
#     mod.mu = pyo.Param(initialize = mu)

#     def rg_stage_cost_rule(m, i):
#         if i == fe_set.at(-1):
#             return pyo.Constraint.Skip
#         else:
#             return (-D_coff*mod.D[i] +
#                     Qr_coff*mod.Qr[i+h]# -
#                     # mod.mu*(pyo.log(mod.x[i, 42]-0.98))
#                     )
#     mod.rg_econ_stage_cost = pyo.Expression(fe_set, rule=rg_stage_cost_rule)


# x_coff = 5.0E1
# D_coff = 1.0E0
# Qr_coff = 1.0E0

# def state_cons_bc_stage_cost(mod):
#     print("No constraint in this stage cost")

# def dist_econ_stage_cost(mod, fe_set, dynamic,):
#     h = mod.h.value

#     def econ_stage_cost_rule(m, i):
#         if dynamic:
#             if i == fe_set.at(-1):
#                 return pyo.Constraint.Skip
#             input_index = i+h
#         else:
#             input_index = i

#         return (-x_coff*m.x[i, 42] - D_coff*m.D[i] + Qr_coff*m.Qr[input_index])

#     mod.econ_stage_cost_expr = pyo.Expression(fe_set,
#                                               rule=econ_stage_cost_rule)


# def construct_stage_cost(mod, fe_set, dynamic):
#     state_cons_bc_stage_cost(mod)
#     dist_econ_stage_cost(mod, fe_set, dynamic)


# def construct_rg_stage_cost(mod, fe_set):
#     h = mod.h.value
#     # mod.mu = pyo.Param(initialize = mu)

#     def rg_stage_cost_rule(m, i):
#         if i == fe_set.at(-1):
#             return pyo.Constraint.Skip
#         else:
#             return (-x_coff*m.x[i, 42] - D_coff*m.D[i] + Qr_coff*m.Qr[i+h]
#                     )
#     mod.rg_econ_stage_cost = pyo.Expression(fe_set, rule=rg_stage_cost_rule)
