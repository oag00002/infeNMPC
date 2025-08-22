from pyomo.environ import (
    Set,
    Var,
    Param,
    Expression,
    Constraint,
    SolverFactory,
#     Objective,
    ComponentUID,
    Reference,
    Block
)
from pyomo.common.collections import ComponentSet, ComponentMap
from pyomo.core.expr.sympy_tools import (
    PyomoSympyBimap,
    Pyomo2SympyVisitor,
    Sympy2PyomoVisitor,
    )
from workspace.eNMPC_red_grad.diehl_dist_model_rules_rg import (
    construct_rg_model_expr
    )
from workspace.models.Diehl_dist.diehl_dist_econ_stage_cost_bound import (
    calibrated_bound_expression_for_KKT_conditions,
)
from workspace.solve_NLP import(
    solve_problem,
)

from workspace.selfStabilizing_eNMPC_utils import (
    create_gradient_matrix_w_pyomo_expression_v2,
    fill_vector_w_pyo_var_at_tp,
    construct_PyomoSympyBimap_with_old_at_diff_var_ind,
)
import numpy as np
import time
import json
from os.path import join


str_states = ("x", "M")
str_controls = ("Rec", "Qr")

def construct_KKT_expressions(olnmpc):
    olnmpc.rg_block = Block()
    b = olnmpc.rg_block
    b.deactivate()
    construct_rg_model_expr(b)

    # Categorize --------------------------------------------------
    # List of differential variables
    DVars = []
    for j in olnmpc.tray:
        for str_state in str_states:
            var = olnmpc.find_component(str_state)
            DVars.append(Reference(var[:,j]))

    # List of algebric variables
    str_AVars = ["dummy_x", "T", "Tdot", "pm", "pn",
                 "V", "L", "y", "hl", "hv",
                 "Vm", "Mv",]
    AVars = []
    for tr in olnmpc.tray:
        for avar in str_AVars:
            if tr == 1 and avar == "Mv":
                AVars.append(Reference(olnmpc.Mv1[:]))
            elif tr == 42 and avar in ["V", "y", "hv"]:
                pass
            elif tr == 42 and avar == "Mv":
                AVars.append(Reference(olnmpc.Mvn[:]))
            else:
                var = olnmpc.find_component(avar)
                AVars.append(Reference(var[:, tr]))
    AVars.append(Reference(olnmpc.D[:]))
    AVars.append(Reference(olnmpc.Qc[:]))

    # List of differential equations
    str_DEqus = ["rg_de_x", "rg_de_M"]
    DEqus = []
    for tr in olnmpc.tray:
        for dequ in str_DEqus:
            equ = b.find_component(dequ)
            indexequ = Reference(equ[:, tr])
            DEqus.append(indexequ)

    # List of algebric equations
    str_AEqus = ["rg_for_dummy_x", "rg_dp", "rg_lTdot", "rg_lpself", "rg_lpn",
                 "rg_gh", "rg_hyd", "rg_gy", "rg_hkl", "rg_hkv",
                 "rg_dvself", "rg_dMV"]
    AEqus = []
    for tr in olnmpc.tray:
        for aequ in str_AEqus:
            if tr == 1 and aequ == "rg_gy":
                AEqus.append(Reference(b.rg_gy0[:]))
            elif tr == 42 and aequ in ["rg_gy", "rg_hkv"]:
                pass
            else:
                equ = b.find_component(aequ)
                AEqus.append(Reference(equ[:,tr]))
    AEqus.append(Reference(b.rg_hrc[:]))

    Phi = [olnmpc.econ_stage_cost]
    controls = [olnmpc.find_component(str_control) for str_control in str_controls]

    nx = len(DVars)
    ny = len(AVars)
    nu = len(controls)
    olnmpc.nx = Param(initialize = nx)
    olnmpc.ny = Param(initialize = ny)
    olnmpc.nu = Param(initialize = nu)
    olnmpc.x_ind = Set(initialize = [i for i in range(olnmpc.nx.value)])
    olnmpc.y_ind = Set(initialize = [i for i in range(olnmpc.ny.value)])
    olnmpc.u_ind = Set(initialize = [i for i in range(olnmpc.nu.value)])

    # Construct KKT conditions at t0 ------------------------------------
    t0 = olnmpc.t.first()
    DVars_t0 = [var[t0] for var in DVars]
    AVars_t0 = [var[t0] for var in AVars]
    DEqus_t0 = [expr[t0] for expr in DEqus]
    AEqus_t0 = [expr[t0] for expr in AEqus]
    Phi_t0 = [Phi[0][t0]]
    controls_t0 = [var[olnmpc.h] for var in controls] # The first control is at t=h

    # Gradient at t0
    dfdx_t0 = create_gradient_matrix_w_pyomo_expression_v2(
        DEqus_t0,
        DVars_t0,
        method="symbolic")
    dfdy_t0 = create_gradient_matrix_w_pyomo_expression_v2(
        DEqus_t0,
        AVars_t0,
        method="symbolic")
    dfdu_t0 = create_gradient_matrix_w_pyomo_expression_v2(
        DEqus_t0,
        controls_t0,
        method="symbolic")

    dgdx_t0 = create_gradient_matrix_w_pyomo_expression_v2(
        AEqus_t0,
        DVars_t0,
        method="symbolic")
    dgdy_t0 = create_gradient_matrix_w_pyomo_expression_v2(
        AEqus_t0,
        AVars_t0,
        method="symbolic")
    dgdu_t0 = create_gradient_matrix_w_pyomo_expression_v2(
        AEqus_t0,
        controls_t0,
        method="symbolic")

    dphidx_t0 = create_gradient_matrix_w_pyomo_expression_v2(
        Phi_t0,
        DVars_t0,
        method="symbolic")
    dphidy_t0 = create_gradient_matrix_w_pyomo_expression_v2(
        Phi_t0,
        AVars_t0,
        method="symbolic")
    dphidu_t0 = create_gradient_matrix_w_pyomo_expression_v2(
        Phi_t0,
        controls_t0,
        method="symbolic")

    nh = 2
    olnmpc.nh = Param(initialize=nh)
    olnmpc.h_ind = Set(initialize=[i for i in range(nh)])

    mu_x42, b.h_x42, mu_x1, b.h_x1 = \
        calibrated_bound_expression_for_KKT_conditions(olnmpc)

    h_t0 = [b.h_x42[0], b.h_x1[0]]
    dhdx_t0 = create_gradient_matrix_w_pyomo_expression_v2(
        h_t0,
        DVars_t0,
        method="symbolic")
    dhdy_t0 = create_gradient_matrix_w_pyomo_expression_v2(
        h_t0,
        AVars_t0,
        method="symbolic")
    dhdu_t0 = create_gradient_matrix_w_pyomo_expression_v2(
        h_t0,
        controls_t0,
        method="symbolic")

    olnmpc.lamb = Var(olnmpc.x_ind, initialize = 10.0)
    olnmpc.nue = Var(olnmpc.y_ind, initialize = 20.0)
    olnmpc.eta = Var(olnmpc.h_ind, initialize = 30.0, bounds=(0.0, None))

    olnmpc.KKT1_var = Var(olnmpc.Nfe, olnmpc.x_ind, initialize = 40.0)
    olnmpc.KKT2_var = Var(olnmpc.Nfe, olnmpc.y_ind, initialize = 50.0)
    olnmpc.KKT3_var = Var(olnmpc.Nfe, olnmpc.u_ind, initialize = 60.0)
    olnmpc.KKT4_var = Var(olnmpc.Nfe, olnmpc.h_ind, initialize = 70.0)

    vec_lamb = np.zeros((nx, 1), dtype = object)
    vec_nue = np.zeros((ny, 1), dtype = object)
    vec_eta = np.zeros((nh, 1), dtype = object)

    vec_KKT1_var_t0 = np.zeros((nx, 1), dtype = object)
    vec_KKT2_var_t0 = np.zeros((ny, 1), dtype = object)
    vec_KKT3_var_t0 = np.zeros((nu, 1), dtype = object)
    vec_KKT4_var_t0 = np.zeros((nh, 1), dtype = object)

    fill_vector_w_pyo_var_at_tp(None, olnmpc.lamb, vec_lamb)
    fill_vector_w_pyo_var_at_tp(None, olnmpc.nue, vec_nue)
    fill_vector_w_pyo_var_at_tp(None, olnmpc.eta, vec_eta)

    fill_vector_w_pyo_var_at_tp(t0, olnmpc.KKT1_var, vec_KKT1_var_t0)
    fill_vector_w_pyo_var_at_tp(t0, olnmpc.KKT2_var, vec_KKT2_var_t0)
    fill_vector_w_pyo_var_at_tp(t0, olnmpc.KKT3_var, vec_KKT3_var_t0)
    fill_vector_w_pyo_var_at_tp(t0, olnmpc.KKT4_var, vec_KKT4_var_t0)

    # KKT conditions for the ss problem
    KKT1_LHS = dphidx_t0 + dfdx_t0.dot(vec_lamb) + dgdx_t0.dot(vec_nue) + dhdx_t0.dot(vec_eta)
    def KKT1_t0_rule(m,i):
        return vec_KKT1_var_t0[i][0] == KKT1_LHS[i][0]
    olnmpc.KKT1_t0 = Constraint(olnmpc.x_ind, rule = KKT1_t0_rule)

    KKT2_LHS = dphidy_t0 + dfdy_t0.dot(vec_lamb) + dgdy_t0.dot(vec_nue) + dhdy_t0.dot(vec_eta)
    def KKT2_t0_rule(m,i):
        return vec_KKT2_var_t0[i][0] == KKT2_LHS[i][0]
    olnmpc.KKT2_t0 = Constraint(olnmpc.y_ind, rule = KKT2_t0_rule)

    KKT3_LHS = dphidu_t0 + dfdu_t0.dot(vec_lamb) + dgdu_t0.dot(vec_nue) + dhdu_t0.dot(vec_eta)
    def KKT3_t0_rule(m,i):
        return vec_KKT3_var_t0[i][0] == KKT3_LHS[i][0] # This is an assignment
    olnmpc.KKT3_t0 = Constraint(olnmpc.u_ind, rule = KKT3_t0_rule)

    negH = np.array([[-b.h_x42[0], 0.],
                     [0., -b.h_x1[0]]])
    # KKT4_LHS = negH.dot(vec_eta_t0)
    mu_list = [mu_x42, mu_x1]
    def KKT4_t0_rule(m,i):
        # return KKT4_LHS[i][0] - olnmpc.mu == 0.0
        return vec_KKT4_var_t0[i][0] == m.eta[i]*negH[i][i] - mu_list[i]
    olnmpc.KKT4_t0 = Constraint(olnmpc.h_ind, rule = KKT4_t0_rule)

    # time1 = time.time()

def initialize_KKT_variables_with_ipopt(olnmpc, csolver):
    # Initialize lamb, nue, eta
    for var in olnmpc.component_objects(Var):
        try:
            var.fix()
        except AttributeError:
            var[...].fix()

    for con in olnmpc.component_objects(Constraint):
        con.deactivate()

    t0 = olnmpc.t.first()
    olnmpc.KKT1_var[t0,:].fix(0)
    olnmpc.KKT2_var[t0,:].fix(0)
    olnmpc.KKT4_var[t0,:].fix(0)
    olnmpc.lamb.unfix()
    olnmpc.nue.unfix()
    olnmpc.eta.unfix()
    olnmpc.KKT1_t0.activate()
    olnmpc.KKT2_t0.activate()
    olnmpc.KKT4_t0.activate()

    solve_problem("initialize lambda, nue, and eta", csolver, olnmpc)

    # time2 = time.time()

    olnmpc.KKT3_var.unfix()
    olnmpc.KKT3_t0.activate()

    solve_problem("initialize KKT3_var", csolver, olnmpc)

    # time3 = time.time()

    olnmpc.x[0, 42].pprint()
    olnmpc.Qr[0].pprint()
    olnmpc.econ_stage_cost[0].pprint()
    Reference(olnmpc.KKT3_var[0, :]).pprint()

    for var in olnmpc.component_objects(Var):
        try:
            var.unfix()
        except AttributeError:
            var[...].unfix()

    for con in olnmpc.component_objects(Constraint):
        con.activate()

    # Initialize KKT1_var, KKT2_var, KKT3_var, KKT4_var @ other points
    var_set_map = ComponentMap([
        (olnmpc.KKT1_var, olnmpc.x_ind),
        (olnmpc.KKT2_var, olnmpc.y_ind),
        (olnmpc.KKT3_var, olnmpc.u_ind),
        (olnmpc.KKT4_var, olnmpc.h_ind)
        ])
    for var, indset in var_set_map.items():
        t0 = olnmpc.Nfe.first()
        for ind in indset:
            tval = var[(t0, ind)].value
            for tind in olnmpc.Nfe:
                if tind != t0:
                    var[(tind, ind)].set_value(tval)

def construct_KKT_constraints_for_stages(olnmpc):
    # Go through the olnmpc.Nfe for each KKT's each non-time indices
    olnmpc.KKT1 = Constraint(olnmpc.Nfe, olnmpc.x_ind)
    olnmpc.KKT2 = Constraint(olnmpc.Nfe, olnmpc.y_ind)
    olnmpc.KKT3 = Constraint(olnmpc.Nfe, olnmpc.u_ind)
    olnmpc.KKT4 = Constraint(olnmpc.Nfe, olnmpc.h_ind)
    KKT_t0_gen_map = ComponentMap([(olnmpc.KKT1_t0, olnmpc.KKT1),
                                   (olnmpc.KKT2_t0, olnmpc.KKT2),
                                   (olnmpc.KKT3_t0, olnmpc.KKT3),
                                   (olnmpc.KKT4_t0, olnmpc.KKT4),
                                   ])

    controls = [olnmpc.find_component(str_control) for str_control in str_controls]
    for con_t0, con_gen in KKT_t0_gen_map.items():
        scalar_set = con_t0.index_set()
        for ind in scalar_set:
            src_body = con_t0[ind].body

            object_map = PyomoSympyBimap()
            visitor = Pyomo2SympyVisitor(object_map)
            sympy_expr = visitor.walk_expression(src_body)

            for tind in olnmpc.Nfe:
                if tind != olnmpc.Nfe.at(-1):
                    index_to_change = ComponentMap([(olnmpc.t, tind), (olnmpc.Nfe, tind)])
                    new_PyomoSympyBimap = construct_PyomoSympyBimap_with_old_at_diff_var_ind(
                                                object_map,
                                                index_in_set_to_change=index_to_change,
                                                control_vars=controls,
                                                sample_time=olnmpc.h.value,
                                                                                            )
                    new_visitor = Sympy2PyomoVisitor(new_PyomoSympyBimap)
                    new_body = new_visitor.walk_expression(sympy_expr)

                    con_gen.add(index=(tind, ind),
                                expr=new_body==0.0)

    del olnmpc.KKT1_t0
    del olnmpc.KKT2_t0
    del olnmpc.KKT3_t0
    del olnmpc.KKT4_t0
    # olnmpc.KKT1_t0.deactivate()
    # olnmpc.KKT2_t0.deactivate()
    # olnmpc.KKT3_t0.deactivate()
    # olnmpc.KKT4_t0.deactivate()

    # time4 = time.time()
    # print("")
    # print("algebraic calculations: ", time1-time0)
    # print("initialize lambda and nue: ", time2-time1)
    # print("initialize lambda, nu, and KKT3_var: ", time3-time2)
    # print("create constraints: ", time4-time3)

def initialize_dual_variables_by_file(olnmpc):
    # fileDir = "/home/kuanhanl/KHL/Research/eNMPC/eNMPC_reduced_gradient/Distillation_column/Check_KKT_conditions"
    fileDir = "Check_KKT_conditions"
    fileName = "lamb_nue_eta_info1.json"
    filePath = join(fileDir, fileName)
    with open(filePath, "r") as fr:
        dual_info = json.load(fr)
    for key, val in dual_info.items():
        cuid = ComponentUID(key)
        var = olnmpc.find_component(cuid)
        var.set_value(val)
