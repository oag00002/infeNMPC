from pyomo.core.expr.calculus.derivatives import Modes, differentiate
from pyomo.core.expr.sympy_tools import (
    PyomoSympyBimap,
    Pyomo2SympyVisitor,
    Sympy2PyomoVisitor,
    )
import numpy as np

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


def fill_vector_w_pyo_var_at_tp(tp, var, vec):
    nrow, ncol = vec.shape
    assert ncol == 1
    if tp is not None:
        assert len([ele for ele in var[tp, :]]) == nrow
        for ind, ele in enumerate(var[tp, :]):
            vec[ind][0] = ele
    else:
        assert len([ele for ele in var[:]]) == nrow
        for ind, ele in enumerate(var[:]):
            vec[ind][0] = ele


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
