"""
eNMPC CSTR model — isothermal CSTR with consecutive reaction A → B.

States: Ca (mol/L), Cb (mol/L)
MV:     Fa0 — feed rate of A (kmol/h), bounds [10, 20]
Objective: maximize B production via economic stage cost -Fa0*(2*Cb - 0.5)
"""
import pyomo.environ as pyo
from pyomo.dae import DerivativeVar, ContinuousSet
from pyomo.core.expr import exp




def variables_initialize(m):
    """
    Declare all Pyomo variables, sets, and parameters for the system model.

    Adds MV/CV index sets, differential state variables, manipulated inputs,
    derivative variables, setpoints, and initial values to the model.

    Parameters
    ----------
    m : pyo.ConcreteModel
        The Pyomo model to which components are added.

    Returns
    -------
    pyo.ConcreteModel
        The model with all variables and parameters declared.
    """
    # ---- Variable Type Indexing ----
    m.MV_index = pyo.Set(initialize=["Fa0"])    # Manipulated variables
    m.CV_index = pyo.Set(initialize=["Ca", "Cb"])      # Controlled variables

    # ---- Plot Display Names in markdown format ----
    m.time_display_name = ["Time (h)"]
    m.MV_display_names = ["F_{A0}\;(kmol/h)"]
    m.CV_display_names = ["C_a\;(mol/L)", "C_b\;(mol/L)"]

    # ---- Differential State Variables ----
    m.Ca = pyo.Var(m.time, initialize=0.5, domain=pyo.NonNegativeReals)  # A concentration
    m.Cb = pyo.Var(m.time, initialize=0.5, domain=pyo.NonNegativeReals)  # B concentration

    # ---- Manipulated Inputs ----
    m.Fa0 = pyo.Var(m.time, initialize=12, bounds=(10,20))    # A feed rate

    # ---- Derivative Variables ----
    m.dCadt = DerivativeVar(m.Ca, initialize=0, wrt=m.time)
    m.dCbdt = DerivativeVar(m.Cb, initialize=0, wrt=m.time)

    # ---- Setpoints for CVs ----
    m.setpoints = pyo.Param(m.CV_index, initialize={"Ca": 0.5, "Cb": 0.5})

    # ---- Initial Values for MVs ----
    m.initial_values = pyo.Param(m.CV_index, initialize={"Ca": 0, "Cb": 0})

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
    # ---- Constants ----
    V = 10
    k = 1.2
    caf = 1

    def Ca_balance_rule(m, t):
        if t not in m.dCadt:
            return pyo.Constraint.Skip
        return m.dCadt[t] == m.Fa0[t] / V * (caf - m.Ca[t]) - k*m.Ca[t]
    def Cb_balance_rule(m, t):
        if t not in m.dCbdt:
            return pyo.Constraint.Skip
        return m.dCbdt[t] == m.Fa0[t] / V * (0 - m.Cb[t]) + k*m.Ca[t]

    m.Ca_balance = pyo.Constraint(m.time, rule=Ca_balance_rule)
    m.Cb_balance = pyo.Constraint(m.time, rule=Cb_balance_rule)

    return m


def custom_objective(m, options):
    stage_cost = lambda m, t: -m.Fa0[t]*(2*m.Cb[t] - 0.5)
    return stage_cost


