import pyomo.environ as pyo
from pyomo.dae import DerivativeVar, ContinuousSet
from pyomo.core.expr import exp


def variables_initialize(m):
    # ---- Time Units ----
    m.time_display_name = ["Time (h)"]

    # ---- Variable Type Indexing ----
    m.MV_index = pyo.Set(initialize=["Fa0", "mc"])    # Manipulated variables
    m.MV_display_names = ["F_{A0}\;(kmol/h)", "\dot{m}_c\;(kmol/h)"]
    m.CV_index = pyo.Set(initialize=["Cc", "T"])      # Controlled variables
    m.CV_display_names = ["C_c\;(mol/L)", "T\;(K)"]
    m.DV_index = pyo.Set(initialize=["Fb0", "UA"])    # Disturbance variables

    # ---- Differential State Variables ----
    m.Ca = pyo.Var(m.time, initialize=1.5, domain=pyo.NonNegativeReals)  # A concentration
    m.Cb = pyo.Var(m.time, initialize=1.5, domain=pyo.NonNegativeReals)  # B concentration
    m.Cc = pyo.Var(m.time, initialize=1.5, domain=pyo.NonNegativeReals)  # C concentration
    m.Cm = pyo.Var(m.time, initialize=1.5, domain=pyo.NonNegativeReals)  # Inert/misc concentration
    m.T = pyo.Var(m.time, initialize=297, domain=pyo.NonNegativeReals)   # Reactor temperature

    # ---- Parameters (Disturbances) ----
    m.Fb0 = pyo.Param(initialize=453.6, mutable=False)  # B feed rate
    m.UA = pyo.Param(initialize=7262, mutable=False)    # Heat transfer coefficient * area

    # ---- Manipulated Inputs ----
    m.Fa0 = pyo.Var(m.time, initialize=35, bounds=(10,100))    # A feed rate
    m.mc = pyo.Var(m.time, initialize=450, bounds=(250,1000))  # Coolant flow rate

    # ---- Derivative Variables ----
    m.dCadt = DerivativeVar(m.Ca, wrt=m.time)
    m.dCbdt = DerivativeVar(m.Cb, wrt=m.time)
    m.dCcdt = DerivativeVar(m.Cc, wrt=m.time)
    m.dCmdt = DerivativeVar(m.Cm, wrt=m.time)
    m.dTdt = DerivativeVar(m.T, wrt=m.time)

    # ---- Setpoints for CVs ----
    m.setpoints = pyo.Param(m.CV_index, initialize={"Cc": 5.18, "T": 396.1})

    # ---- Initial Values for MVs ----
    m.initial_values = pyo.Param(m.MV_index, initialize={"Fa0": 36.3, "mc": 453.6})

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
    Ta1 = 288.7     # Ambient coolant inlet temperature [K]
    CpW = 18        # Specific heat capacity of coolant water [cal/mol-K]
    T0 = 297        # Reference temperature [K]
    dH = -20013     # Heat of reaction [cal/mol]
    V = 5           # Reactor volume [L]
    Fm0 = 45.4      # Feed rate of inert/miscellaneous stream [mol/min]

    # ---- Algebraic Equations ----
    m.k = {t: 1.696e13 * exp(-18012 / 1.987 / (m.T[t])) for t in m.time}
    m.ra = {t: 0 - (m.k[t] * m.Ca[t]) for t in m.time}
    m.rb = {t: 0 - (m.k[t] * m.Ca[t]) for t in m.time}
    m.rc = {t: m.k[t] * m.Ca[t] for t in m.time}
    m.Na = {t: m.Ca[t] * V for t in m.time}
    m.Nb = {t: m.Cb[t] * V for t in m.time}
    m.Nc = {t: m.Cc[t] * V for t in m.time}
    m.Nm = {t: m.Cm[t] * V for t in m.time}
    m.ThetaCp = {t: 35 + m.Fb0 / m.Fa0[t] * 18 + Fm0 / m.Fa0[t] * 19.5 for t in m.time}
    m.v0 = {t: m.Fa0[t] / 14.8 + m.Fb0 / 55.3 + Fm0 / 24.7 for t in m.time}
    m.Ta2 = {t: m.T[t] - ((m.T[t] - Ta1) * exp(0 - (m.UA / (CpW * m.mc[t])))) for t in m.time}
    m.Ca0 = {t: m.Fa0[t] / m.v0[t] for t in m.time}
    m.Cb0 = {t: m.Fb0 / m.v0[t] for t in m.time}
    m.Cm0 = {t: Fm0 / m.v0[t] for t in m.time}
    m.Qr2 = {t: m.mc[t] * CpW * (m.Ta2[t] - Ta1) for t in m.time}
    m.Qr1 = {t: m.Fa0[t] * m.ThetaCp[t] * (m.T[t] - T0) for t in m.time}
    m.Qr = {t: m.Qr1[t] + m.Qr2[t] for t in m.time}
    m.Qg = {t: m.ra[t] * V * dH for t in m.time}
    m.tau_tc = {t: V / m.v0[t] for t in m.time}
    m.NCp = {t: m.Na[t] * 35 + m.Nb[t] * 18 + m.Nc[t] * 46 + m.Nm[t] * 19.5 for t in m.time}

    def Ca_balance_rule(m, t):
        return m.dCadt[t] == (1 / m.tau_tc[t]) * (m.Ca0[t] - m.Ca[t]) + m.ra[t]
    def Cb_balance_rule(m, t):
        return m.dCbdt[t] == (1 / m.tau_tc[t]) * (m.Cb0[t] - m.Cb[t]) + m.rb[t]
    def Cc_balance_rule(m, t):
        return m.dCcdt[t] == (1 / m.tau_tc[t]) * (0 - m.Cc[t]) + m.rc[t]
    def Cm_balance_rule(m, t):
        return m.dCmdt[t] == (1 / m.tau_tc[t]) * (m.Cm0[t] - m.Cm[t])
    def energy_balance_rule(m, t):
        return m.dTdt[t] == (m.Qg[t] - m.Qr[t]) / m.NCp[t]

    m.Ca_balance = pyo.Constraint(m.time, rule=Ca_balance_rule)
    m.Cb_balance = pyo.Constraint(m.time, rule=Cb_balance_rule)
    m.Cc_balance = pyo.Constraint(m.time, rule=Cc_balance_rule)
    m.Cm_balance = pyo.Constraint(m.time, rule=Cm_balance_rule)
    m.energy_balance = pyo.Constraint(m.time, rule=energy_balance_rule)

    return m


# ---- Testing the Model ----
if __name__ == '__main__':
    horizon = 10
    m = pyo.ConcreteModel()
    m.time = ContinuousSet(bounds=(0, 10))

    print('Testing Model Equations')
    try:
        m = equations_write(m)
        print('Equation Writing Successful')
    except Exception as e:
        print(f'Equation Writing Failed: {e}')
    model_display_flag = True
    if model_display_flag:
        m.pprint()

    assert isinstance(m, pyo.ConcreteModel)