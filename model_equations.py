import pyomo.environ as pyo
from pyomo.dae import DerivativeVar, ContinuousSet
from pyomo.core.expr import exp
import numpy as np


def variables_initialize(m):
    # ---- Time Units ----
    m.time_display_name = ["Time (s)"]

    # ---- Variable Type Indexing ----
    m.MV_index = pyo.Set(initialize=["F"])    # Manipulated variables
    m.MV_display_names = ["F\;(N)"]
    m.CV_index = pyo.Set(initialize=["theta", "x"])      # Controlled variables
    m.CV_display_names = ["\theta", "x\;(m)"]

    # ---- Differential State Variables ----
    m.x = pyo.Var(m.time, initialize=0, bounds=(-1000, 1000))
    m.x_dot = pyo.Var(m.time, initialize=0, bounds=(-100, 100))
    m.theta = pyo.Var(m.time, initialize=np.radians(10), bounds=(np.radians(-180), np.radians(180)))
    m.theta_dot = pyo.Var(m.time, initialize=0, bounds=(np.radians(-90), np.radians(90)))

    # ---- Parameters (Disturbances) ----
    m.b = pyo.Param(initialize=0.1, mutable=True)

    # ---- Manipulated Inputs ----
    m.F = pyo.Var(m.time, initialize=0, bounds=(-5,5))

    # ---- Derivative Variables ----
    m.dx = DerivativeVar(m.x, wrt=m.time, initialize=0)
    m.dx_dot = DerivativeVar(m.x_dot, wrt=m.time, initialize=0)
    m.dtheta = DerivativeVar(m.theta, wrt=m.time, initialize=0)
    m.dtheta_dot = DerivativeVar(m.theta_dot, wrt=m.time, initialize=0)

    # ---- Setpoints for CVs ----
    m.setpoints = pyo.Param(m.CV_index, initialize={"theta": 0, "x": 0})

    # ---- Initial Values for MVs ----
    m.initial_values = pyo.Param(m.CV_index, initialize={"theta": np.radians(10), "x": 0})

    return m


def equations_write(m):
    """
    Define and apply all variables, differential equations, algebraic expressions,
    and constraints required to generate a model.

    Parameters:
    -----------
    m : pyomo.ConcreteModel
        A Pyomo model object to which variables and constraints are added.

    Returns:
    --------
    m : pyomo.ConcreteModel
        The model with all relevant variables and constraints defined.
    """
    m_p = 0.2
    l = 0.3
    m_c = 0.5
    k = 1 / 3
    g = 9.8

    def x_dyn_rule(m, t):
        return m.dx[t] == m.x_dot[t]
    def x_dot_dyn_rule(m, t):
        return m.dx_dot[t] * (m_p * pyo.cos(m.theta[t])**2 - (1 + k) * (m_c + m_p)) == (m_p * g * pyo.sin(m.theta[t]) * pyo.cos(m.theta[t]) - (1 + k) * (m.F[t] + m_p * l * m.theta_dot[t]**2 * pyo.sin(m.theta[t]) - m.b * m.x_dot[t]))
    def theta_dyn_rule(m, t):
        return m.dtheta[t] == m.theta_dot[t]
    def theta_dot_dyn_rule(m, t):
        return m.dtheta_dot[t] * ((1 + k) * (m_c + m_p) * l - m_p * l * pyo.cos(m.theta[t])**2) == (m_c * g * pyo.sin(m.theta[t]) - pyo.cos(m.theta[t]) * (m.F[t] + m_p * l * m.theta_dot[t]**2 * pyo.sin(m.theta[t])))

    m.x_dyn = pyo.Constraint(m.time, rule=x_dyn_rule)
    m.x_dot_dyn = pyo.Constraint(m.time, rule=x_dot_dyn_rule)
    m.theta_dyn = pyo.Constraint(m.time, rule=theta_dyn_rule)
    m.theta_dot_dyn = pyo.Constraint(m.time, rule=theta_dot_dyn_rule)

    return m


# ---- Testing the Model ----
if __name__ == '__main__':
    horizon = 10
    m = pyo.ConcreteModel()
    m.time = ContinuousSet(bounds=(0, 10))

    print('Testing Model Equations')
    try:
        m = variables_initialize(m)
        m = equations_write(m)
        print('Equation Writing Successful')
    except Exception as e:
        print(f'Equation Writing Failed: {e}')
    model_display_flag = True
    if model_display_flag:
        m.pprint()

    assert isinstance(m, pyo.ConcreteModel)