import re
from pyomo.dae import DerivativeVar


def _get_variable_key_for_data(model, name):
    """
    Generate a string key like 'x[1,2,*]' from a variable name string.

    Parameters
    ----------
    model : ConcreteModel
        The Pyomo model containing the variable.
    name : str
        Variable name with optional index, e.g., 'x[1,2]'.

    Returns
    -------
    str
        A string key formatted for get_data_from_key, e.g., 'x[*]' or 'x[1,2,*]'.
    """
    var_name, index = _parse_indexed_name(name)
    var_obj = getattr(model, var_name)
    base_name = var_obj.name.split("[")[0]
    if not index:
        return f"{base_name}[*]"
    index_str = ",".join(str(i) for i in index) + ",*"
    return f"{base_name}[{index_str}]"


def _parse_indexed_name(name):
    """
    Given 'x[1,2]', returns ('x', (1,2)).
    If 'Qr', returns ('Qr', ()).
    If 'Qr[*]', returns ('Qr', ('*',)).
    """
    match = re.match(r"^([a-zA-Z_]\w*)(?:\[(.*)\])?$", name)
    if not match:
        raise ValueError(f"Invalid variable format: {name}")
    var_name = match.group(1)
    index_str = match.group(2)
    if index_str:
        index = tuple(i.strip() if i.strip() == "*" else eval(i.strip()) for i in index_str.split(','))
    else:
        index = ()
    return var_name, index


def _add_time_indexed_expression(model, var_name, t):
    """
    Return an expression m.var[i1, ..., t] for given var_name = 'var[i1, ...]'.

    Parameters
    ----------
    model : ConcreteModel
    var_name : str
    t : time index

    Returns
    -------
    pyomo expression
    """
    name, base_index = _parse_indexed_name(var_name)
    var_obj = getattr(model, name)
    if base_index:
        return var_obj[base_index + (t,)]
    else:
        return var_obj[t]


def _get_disc_eq_time_points(m):
    """Return a set of all unique collocation time points used in discretization equations."""
    time_points = set()
    for var in getattr(m, 'deriv_vars', []):
        var_name = var.getname().split('.')[-1]  # Extract variable name from full path
        disc_eq = getattr(m, f"{var_name}_disc_eq", None)
        if disc_eq is not None:
            for idx in disc_eq:
                time = idx[-1] if isinstance(idx, tuple) else idx
                time_points.add(time)
    return sorted(time_points)


def _get_derivative_and_state_vars(model):
    """
    Return sets of DerivativeVar and state variable components.
    """
    deriv_vars = set()
    state_vars = set()

    for deriv in model.component_objects(DerivativeVar, descend_into=True):
        deriv_vars.add(deriv)
        state_vars.add(deriv.get_state_var())

    return deriv_vars, state_vars