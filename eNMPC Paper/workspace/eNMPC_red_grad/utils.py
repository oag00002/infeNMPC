from pyomo.environ import (
    Expression,
    Constraint,
    TerminationCondition,
    )
from pyomo.environ import value as pyo_val
from pyomo.dae import DerivativeVar
from pyomo.core.expr.calculus.derivatives import Modes, differentiate
# from workspace.utils import organize_pyomo_expr
import numpy as np
from pyomo.common.collections import ComponentSet, ComponentMap
from pyomo.core.expr.sympy_tools import (
    PyomoSympyBimap,
    Pyomo2SympyVisitor,
    Sympy2PyomoVisitor,
    )

def solve_problem(step, solver, m):
    print("")
    print("iteration = ", step, ", Solve " + m.name)
    results = solver.solve(m, tee = True)
    #check results
    assert results.solver.termination_condition == TerminationCondition.optimal


def categorize_flatten_dae_expr_into_diff_alge_lists(flatten_dae,
                                                     diff_set,
                                                     alge_set = None):
    diff_list = list()
    alge_list = list()
    for comp in flatten_dae:
        t0 = comp.index_set().first()
        parent_comp = comp[t0].parent_component()
        if parent_comp in diff_set:
            diff_list.append(comp)
        else:
            if alge_set is None:
                alge_list.append(comp)
            else:
                if parent_comp in alge_set:
                    alge_list.append(comp)
                else:
                    raise RuntimeError("Cannot categorize ", parent_comp.name)

    return diff_list, alge_list


def categorize_flatten_dae_var_into_diff_alge_lists(time,
                                                    flatten_dae,
                                                    diff_set,
                                                    control_set,
                                                    alge_set = None,
                                                    exclude_set = None,
                                                    ):
    diff_list = list()
    alge_list = list()
    control_list = list()
    deri_list = list()

    t0 = time.first()
    if exclude_set is not None:
        exclude_comps_ref = [comp.referent for comp in exclude_set]

    for comp in flatten_dae:
        parent_comp = comp[t0].parent_component()
        if exclude_set is not None:
            if comp.referent in exclude_comps_ref:
                print(comp.referent, " is excluded!")
                continue

        if type(parent_comp) is DerivativeVar:
            deri_list.append(comp)
        elif parent_comp in diff_set:
            diff_list.append(comp)
        elif parent_comp in control_set:
            control_list.append(comp)
        else:
            if alge_set is None:
                alge_list.append(comp)
            else:
                if parent_comp in alge_set:
                    alge_list.append(comp)
                else:
                    raise RuntimeError("Cannot categorize ", parent_comp.name)

    return diff_list, alge_list, control_list, deri_list


def create_gradient_matrix_w_pyomo_expression(b,
                                              expr_info,
                                              var_info,
                                              gradient_name,
                                              var_is_u = False):
    m = b.parent_block()
    if m is None:
        m = b

    expr_ind_set, expr_list = expr_info
    var_ind_set, var_list = var_info

    h = m.h.value
    def gradient_rule(b, i, j, k):
        if k == m.Nfe.at(-1):
            return Constraint.Skip
        else:
            if var_is_u:
                ind = k+h
            else:
                ind = k

            sym_deriv = differentiate(expr_list[j-1][k], var_list[i-1][ind],
                                      mode = Modes.reverse_symbolic)

            return sym_deriv

    b.add_component(gradient_name, Expression(var_ind_set,
                                              expr_ind_set,
                                              m.Nfe,
                                              rule = gradient_rule))

    gradient_comp = b.find_component(gradient_name)

    num_var = len(var_ind_set)
    num_expr = len(expr_ind_set)
    matrix_dict = dict()
    from pyomo.core.base.expression import _GeneralExpressionData
    for k in m.Nfe:
        if k != m.Nfe.at(-1):
            matrix = np.zeros((num_var, num_expr),
                              dtype=_GeneralExpressionData)
            for i in var_ind_set:
                for j in expr_ind_set:
                    matrix[i-1][j-1] = gradient_comp[(i,j,k)]
            matrix_dict[k] = matrix

    return matrix_dict

def create_gradient_matrix_w_pyomo_expression_v2(
                                                expr_list,
                                                var_list,
                                                method = "symbolic",
                                                 ):
    if method == "numeric":
        mode = Modes.reverse_numeric
        dtype = 'float64'
    elif method == "symbolic":
        mode = Modes.reverse_symbolic
        dtype = object
    else:
        raise RuntimeError(method, "is not suppported!")

    matrix = np.empty(shape = [len(var_list), 0],
                      dtype=dtype)
    for expr_ind, expr_ele in enumerate(expr_list):
        sym_deriv = differentiate(expr_ele,
                                  wrt_list = var_list,
                                  mode = mode)
        matrix = np.append(matrix, np.array([sym_deriv]).T, axis=1)

    return matrix


def organize_expr_in_matrix_dict(matrix_dict):
    for key, val in matrix_dict.items():
        nrow, ncol = val.shape
        for i in range(nrow):
            for j in range(ncol):
                ele = val[i][j]
                new_ele = organize_pyomo_expr(ele)
                val[i][j] = new_ele

    return matrix_dict

def organize_expr_in_matrix(matrix):
    nrow, ncol = matrix.shape
    for i in range(nrow):
        for j in range(ncol):
            ele = matrix[i][j]
            new_ele = organize_pyomo_expr(ele)
            matrix[i][j] = new_ele

    return matrix

def evalue_pyomo_expr_in_matrix(input_mat):
    nrow, ncol = input_mat.shape
    output_matrix = np.zeros((nrow, ncol))
    for i in range(nrow):
        for j in range(ncol):
            output_matrix[i][j] = pyo_val(input_mat[i][j])
    return output_matrix

def compare_two_matrix(mat1, mat2):
    # check dim
    nrow1, ncol1 = mat1.shape
    nrow2, ncol2 = mat2.shape
    if (nrow1 != nrow2):
        print("nrow is different")
    if (ncol1 != ncol2):
        print("ncol is different")

    diff_list = list()
    for i in range(nrow1):
        for j in range(ncol1):
            val1 = mat1[i][j]
            val2 = mat2[i][j]
            if (val1-val2) > 1E-5:
                diff_list.append((i, j))

    return diff_list


def construct_PyomoSympyBimap_with_old_at_diff_var_ind(
        old_PyomoSympyBimap,
        index_in_set_to_change,
        control_vars=None,
        sample_time=None,
                                                       ):
    # index_in_set_to_change = ComponentMap{m.t: 2, m.s: 20}

    # First, make sure index to change is in the set
    for set_to_change, target_ind in index_in_set_to_change.items():
        if target_ind not in set_to_change:
            raise RuntimeError(target_ind, "is not in set ", set_to_change.name)

    if control_vars is not None and sample_time is None:
        raise RuntimeError("Controls are given with sampling time!")

    new_PyomoSympyBimap = PyomoSympyBimap()
    for sympy_comp, pyomo_comp in old_PyomoSympyBimap.sympy2pyomo.items():

        parent_comp = pyomo_comp.parent_component()
        subsets = list(parent_comp.index_set().subsets())

        index_update = pyomo_comp.index()
        for set_to_change, target_ind in index_in_set_to_change.items():
            if set_to_change in subsets:
                if type(index_update) is int or type(index_update) is float:
                    assert len(subsets) == 1
                    if parent_comp in control_vars:
                        index_update = target_ind + sample_time
                    else:
                        index_update = target_ind
                    break
                elif type(index_update) is tuple:
                    ind = subsets.index(set_to_change)
                    temp = list(index_update)
                    if parent_comp in control_vars:
                        temp[ind] = target_ind + sample_time
                    else:
                        temp[ind] = target_ind
                    index_update = tuple(temp)
                else:
                    raise RuntimeError("type of index is neither int nor tuple")

        new_gen_var = parent_comp[index_update]

        new_PyomoSympyBimap.sympy2pyomo[sympy_comp] = new_gen_var
        new_PyomoSympyBimap.pyomo2sympy[new_gen_var] = sympy_comp

    new_PyomoSympyBimap.i = old_PyomoSympyBimap.i

    return new_PyomoSympyBimap


def change_index_in_expression(
        src_expr,
        new_index_info,
        control_vars=None,
        sample_time=None,
                               ):
    object_map = PyomoSympyBimap()
    visitor = Pyomo2SympyVisitor(object_map)
    sympy_expr = visitor.walk_expression(src_expr)

    index_in_set_to_change = ComponentMap(new_index_info)
    new_PyomoSympyBimap = construct_PyomoSympyBimap_with_old_at_diff_var_ind(
                                object_map,
                                index_in_set_to_change=index_in_set_to_change,
                                control_vars=control_vars,
                                sample_time=sample_time,
                                                                    )
    new_visitor = Sympy2PyomoVisitor(new_PyomoSympyBimap)
    new_expr = new_visitor.walk_expression(sympy_expr)

    return new_expr