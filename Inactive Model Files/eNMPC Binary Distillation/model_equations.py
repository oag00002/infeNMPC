import pyomo.environ as pyo
from pyomo.dae import DerivativeVar, ContinuousSet
from pyomo.core.expr import sqrt


def variables_initialize(m):
    print('Initializing Variables')

    # ---- Time Units ----
    m.time_display_name = ["Time (h)"]

    # ---- Variable Type Indexing ----
    m.MV_index = pyo.Set(initialize=["Qr", "Rec"])    # Manipulated variables
    m.MV_display_names = ["Q_r", "Rec"]
    m.CV_index = pyo.Set(initialize=["x[1]", "T[29]"])      # Controlled variables
    m.CV_display_names = ["x_1", "T_{29}"]
    m.DV_index = pyo.Set(initialize=["xf"])    # Disturbance variables


    m.Ntray = pyo.Param(initialize=42)
    m.tray = pyo.Set(initialize=range(1, m.Ntray+1))

    m.feedTray = pyo.Param(initialize=21)
    m.feed = pyo.Param(
        m.tray,
        initialize=lambda m, k: 57.5294 if k == m.feedTray else 0.0,
        mutable=True,
    )


    def _M_init(m, k, i):
        if k == 1:
            return 2*105500.0
        elif 2 <= k <= m.feedTray:
            return 3340.0+(3772.0-3340.0)/(m.feedTray-2)*(k-2)
        else:
            return 2890.0+(4650.0-2890.0)/(m.Ntray-m.feedTray-1)*(k-m.feedTray-1)


    def _x_init(m, k, i):
        if 1 <= k <= 16:
            return 0.999*k/m.Ntray
        else:
            return 0.36+(0.98-0.36)/(m.Ntray-m.feedTray-1)*(k-m.feedTray-1)
        
    
    def _y_init(m, k, i):
        turnpt = 11
        if 1 <= k <= turnpt:
            return 0.064+(0.56-0.064)/(turnpt-1)*(k-1)
        elif turnpt < k <= m.feedTray-1:
            return 0.58+(0.64-0.58)/(m.feedTray-turnpt-2)*(k-turnpt-1)
        else:
            return 0.64+(0.99-0.67)/(m.Ntray-m.feedTray)*(k-m.feedTray)
        
    def _T_init(m, k, i):
        return 336.0 + (370.0 - 336.0)/(m.Ntray-1)*(k-1)

    # ---- Differential State Variables ----
    m.M = pyo.Var(m.tray, m.time, initialize=_M_init)   # Liquid hold-up
    m.x = pyo.Var(m.tray, m.time, initialize=_x_init, bounds=(0, 1))    # Mole-fraction

    # ---- Derivative Variables ----
    m.Mdot = DerivativeVar(m.M, initialize=0.0)
    m.xdot = DerivativeVar(m.x, initialize=0.0)

    m.setpoint_index = m.CV_index.union(m.MV_index)  # Setpoints for controlled and manipulated variables

    # ---- Setpoints for CVs ----
    m.setpoints = pyo.Param(m.setpoint_index, initialize={"x[1]": 0.002170937656550695, "T[29]": 354.3311105693665, "Qr": 1.4410950565834773, "Rec": 0.22732221769393682})

    # ---- Initial Values for MVs ----
    m.initial_values = pyo.Param(m.MV_index, initialize={"Qr": 1.5, "Rec": 10})

    # Tray temperature
    m.T = pyo.Var(m.tray, m.time, initialize=_T_init, bounds=(200, 450))
    m.Tdot = pyo.Var(m.tray, m.time, initialize=1e-05)  # not really a der_var

    # Vapor mole flowrate
    m.V = pyo.Var(m.tray, m.time, initialize=44.0, bounds=(0, 1e3))

    # Vapor mole fraction
    m.y = pyo.Var(m.tray, m.time, initialize=_y_init, bounds=(0, 1))

    # Condenser heat duty
    m.Qc = pyo.Var(m.time, initialize=1.6, bounds=(0, 1e2))
    # Distillate
    m.D = pyo.Var(m.time, initialize=18.33, bounds = (0, 1e4))


    # Reflux ratio
    m.Rec = pyo.Var(m.time, initialize=1.20, bounds=(0.1, 99.999))
    # Re-boiler heat duty
    m.Qr = pyo.Var(m.time, initialize=1.65, bounds=(0.0, None))

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
    print('Writing Model Equations')

    def _p_init(m, k):
        ptray = 9.39e+04
        if k <= m.feedTray-1:
            return _p_init(m, m.feedTray) + m.pstrip*(m.feedTray - k)
        elif m.feedTray-1 < k < m.Ntray:
            return ptray + m.prect*(m.Ntray - k)
        elif k == m.Ntray:
            return ptray
    
    def clipping(x):
        # Smooth approximation of the max(x,0) function
        return 0.5*(x+sqrt(x**2 + 0.001))

    # Feed mole fraction
    m.xf = pyo.Param(initialize=0.32, mutable=True)
    # Feed enthalpy
    m.hf = pyo.Param(initialize=9081.3)

    m.hlm0 = pyo.Param(initialize=2.6786e-04)
    m.hlma = pyo.Param(initialize=-0.14779)
    m.hlmb = pyo.Param(initialize=97.4289)
    m.hlmc = pyo.Param(initialize=-2.1045e04)

    m.hln0 = pyo.Param(initialize=4.0449e-04)
    m.hlna = pyo.Param(initialize=-0.1435)
    m.hlnb = pyo.Param(initialize=121.7981)
    m.hlnc = pyo.Param(initialize=-3.0718e04)

    m.r = pyo.Param(initialize=8.3147)
    m.a = pyo.Param(initialize=6.09648)
    m.b = pyo.Param(initialize=1.28862)
    m.c1 = pyo.Param(initialize=1.016)
    m.d = pyo.Param(initialize=15.6875)
    m.l = pyo.Param(initialize=13.4721)
    m.f = pyo.Param(initialize=2.615)

    m.gm = pyo.Param(initialize=0.557)
    m.Tkm = pyo.Param(initialize=512.6)
    m.Pkm = pyo.Param(initialize=8.096e06)

    m.gn = pyo.Param(initialize=0.612)
    m.Tkn = pyo.Param(initialize=536.7)
    m.Pkn = pyo.Param(initialize=5.166e06)

    m.CapAm = pyo.Param(initialize=23.48)
    m.CapBm = pyo.Param(initialize=3626.6)
    m.CapCm = pyo.Param(initialize=-34.29)

    m.CapAn = pyo.Param(initialize=22.437)
    m.CapBn = pyo.Param(initialize=3166.64)
    m.CapCn = pyo.Param(initialize=-80.15)

    m.pstrip = pyo.Param(initialize=250)
    m.prect = pyo.Param(initialize=190)

    m.p = pyo.Param(m.tray, initialize=_p_init)
    m.alpha = pyo.Param(
        m.tray, initialize=lambda m, k: 0.62 if k <= m.feedTray else 0.35
    )

    # Substitution for Algebraic constraints
    Mv = {} # Liquid Volume Holdup
    Vm = {} # Molar Volume
    pm = {} # Vapor Pressure by Antoine's equation for methanol
    pn = {} # Vapoi Pressure by Antoine's equation for n-propanol
    hl = {} # Liquid Enthalpy
    hv = {} # Vapor Enthalpy
    L = {} # Liquid flowrate
    D = {} # Distillate flowrate

    Mdot_augmented = {}
    Tdot_augmented = {}
    xdot_augmented = {}
    for i in m.time:
        for k in m.tray:
            Vm[k, i] = m.x[k, i]*((1/2288)*0.2685**(1 + (1 - m.T[k, i]/512.4)**0.2453))+(1 - m.x[k, i])*((1/1235)*0.27136**(1 + (1 - m.T[k, i]/536.4)**0.24))
            # Vm[k, i] = clipping(Vm[k, i])
            pm[k, i] = pyo.exp(m.CapAm - m.CapBm / (m.T[k, i] + m.CapCm))
            pn[k, i] = pyo.exp(m.CapAn - m.CapBn / (m.T[k, i] + m.CapCn))
            Mv[k, i] = Vm[k, i] * m.M[k, i]

            hl[k, i] = m.x[k, i]*(m.hlm0*m.T[k, i]**3 + m.hlma*m.T[k, i]**2 + m.hlmb*m.T[k, i] + m.hlmc) + (1 - m.x[k, i])*(m.hln0*m.T[k, i]**3 + m.hlna*m.T[k, i]**2 + m.hlnb*m.T[k, i] + m.hlnc)
            hv[k, i] = m.y[k, i]*(m.hlm0*m.T[k, i]**3 + m.hlma*m.T[k, i]**2 + m.hlmb*m.T[k, i] + m.hlmc + m.r*m.Tkm*pyo.sqrt(1-(m.p[k]/m.Pkm)*(m.Tkm/m.T[k, i])**3)*(m.a-m.b*m.T[k, i]/m.Tkm + m.c1*(m.T[k, i]/m.Tkm)**7 + m.gm*(m.d-m.l*m.T[k, i]/m.Tkm + m.f*(m.T[k, i]/m.Tkm)**7))) + (1 - m.y[k, i])*(m.hln0*m.T[k, i]**3 + m.hlna*m.T[k, i]**2 + m.hlnb*m.T[k, i] + m.hlnc + m.r*m.Tkn*pyo.sqrt(1-(m.p[k]/m.Pkn)*(m.Tkn/m.T[k, i])**3)*(m.a-m.b*m.T[k, i]/m.Tkn + m.c1*(m.T[k, i]/m.Tkn)**7 + m.gn*(m.d - m.l*m.T[k, i]/m.Tkn + m.f*(m.T[k, i]/m.Tkn)**7)))
            if k == 1:
                L[k, i] = 0.166 * (clipping(Mv[k, i] - 8.5) )** 1.5/Vm[k, i]
            elif k == m.Ntray:
                L[k, i] = 0.166 * (clipping(Mv[k, i] - 0.17) )** 1.5/Vm[k, i]
            else:
                L[k, i] = 0.166 * (clipping(Mv[k, i] - 0.155) )** 1.5/Vm[k, i]
            global_time_scale = 1/60  # Convert from minutes to hours
            
            Mdot_augmented[k, i] = m.Mdot[k, i]*global_time_scale
            Tdot_augmented[k, i] = m.Tdot[k, i]*global_time_scale
            xdot_augmented[k, i] = m.xdot[k, i]*global_time_scale

        D[i] = L[m.Ntray + 0, i] / m.Rec[i]

    # Overall mass balances
    def _de_M_rule(m, k, i):
        if k == 1:
            return Mdot_augmented[1, i] == \
                (L[2, i] - L[1, i] - m.V[1, i])
        elif k == m.Ntray:
            return Mdot_augmented[m.Ntray + 0, i] == \
                (m.V[m.Ntray - 1, i] - L[k, i] - D[i])
        else:
            return Mdot_augmented[k, i] == \
                (m.V[k - 1, i] - m.V[k, i] +
                 L[k + 1, i] - L[k, i] +
                 m.feed[k])

    # Component mass balance
    def _de_x_rule(m, k, i):
        if k == 1:
            return xdot_augmented[1, i] * m.M[1, i] == \
                (L[2, i] * (m.x[2, i] - m.x[1, i]) -
                 m.V[1, i] * (m.y[1, i] - m.x[1, i]))
        elif k == m.Ntray:
            return xdot_augmented[m.Ntray + 0, i] * m.M[m.Ntray, i] == \
                (m.V[m.Ntray - 1, i] * (m.y[m.Ntray - 1, i] - m.x[m.Ntray, i]))
        else:
            return xdot_augmented[k, i] * m.M[k, i] == \
                (m.V[k - 1, i] * (m.y[k - 1, i] - m.x[k, i]) +
                 L[k + 1, i] * (m.x[k + 1, i] - m.x[k, i]) -
                 m.V[k, i] * (m.y[k, i] - m.x[k, i]) +
                 m.feed[k] * (m.xf - m.x[k, i]))


    # Energy balance
    def _gh_rule(m, k, i):
        if k == 1:
            return m.M[1, i] * (
                    xdot_augmented[1, i] * (
                    (m.hlm0 - m.hln0) * m.T[1, i] ** 3 +
                    (m.hlma - m.hlna) * m.T[1, i] ** 2 +
                    (m.hlmb - m.hlnb) * m.T[1, i] +
                    (m.hlmc - m.hlnc)
            ) +
                    Tdot_augmented[1, i] * (
                            3 * m.hln0 * m.T[1, i] ** 2 +
                            2 * m.hlna * m.T[1, i] +
                            m.hlnb +
                            m.x[1, i] * (
                                    3 * (m.hlm0 - m.hln0) * m.T[1, i] ** 2 +
                                    2 * (m.hlma - m.hlna) * m.T[1, i] +
                                    (m.hlmb - m.hlnb)
                            )
                    )
            ) == \
                (L[2, i] * (hl[2, i] -  hl[1, i]) -
                 m.V[1, i] * (hv[1, i] - hl[1, i]) +
                 m.Qr[i]*1.0e6
                 )

        elif k == m.Ntray:
            return m.M[m.Ntray, i] * (
                    xdot_augmented[m.Ntray + 0, i] * (
                    (m.hlm0 - m.hln0) * m.T[m.Ntray, i] ** 3 +
                    (m.hlma - m.hlna) * m.T[m.Ntray, i] ** 2 +
                    (m.hlmb - m.hlnb) * m.T[m.Ntray, i] +
                    (m.hlmc - m.hlnc)
            ) +
                    Tdot_augmented[m.Ntray + 0, i] * (
                            3 * m.hln0 * m.T[m.Ntray, i] ** 2 +
                            2 * m.hlna * m.T[m.Ntray, i] +
                            m.hlnb +
                            m.x[m.Ntray, i] * (
                                    3 * (m.hlm0 - m.hln0) * m.T[m.Ntray, i] ** 2 +
                                    2 * (m.hlma - m.hlna) * m.T[m.Ntray, i] +
                                    (m.hlmb - m.hlnb)
                            )
                    )
            ) == \
                (m.V[m.Ntray - 1, i] * (hv[m.Ntray - 1, i] - hl[m.Ntray + 0, i]) -
                 m.Qc[i]*1e06
                 )
        else:
            return m.M[k, i] * (
                    xdot_augmented[k, i] * (
                    (m.hlm0 - m.hln0) * (m.T[k, i] ** 3) +
                    (m.hlma - m.hlna) * (m.T[k, i] ** 2) +
                    (m.hlmb - m.hlnb) * m.T[k, i] +
                    (m.hlmc - m.hlnc)
            ) +
                    Tdot_augmented[k, i] * (
                            3 * m.hln0 * (m.T[k, i] ** 2) +
                            2 * m.hlna * m.T[k, i] +
                            m.hlnb +
                            m.x[k, i] * (
                                    3 * (m.hlm0 - m.hln0) * (m.T[k, i] ** 2) +
                                    2 * (m.hlma - m.hlna) * m.T[k, i] +
                                    (m.hlmb - m.hlnb)
                            )
                    )
            ) == \
                (m.V[k - 1, i] * (hv[k - 1, i] - hl[k, i]) +
                 L[k + 1, i] * (hl[k + 1, i] - hl[k, i]) -
                 m.V[k, i] * (hv[k, i] - hl[k, i]) +
                 m.feed[k] * (m.hf - hl[k, i])
                 )


    # Raoult's law
    def _dp_rule(m, k, i):
        return m.p[k] == m.x[k, i] * pm[k, i] + (1 - m.x[k, i]) * pn[k, i]

    # Derivative of T for index reduction
    def _lTdot_rule(m, k, i):
        return global_time_scale*m.Tdot[k, i] == \
            -(pm[k, i] - pn[k, i]) * global_time_scale*m.xdot[k, i] / \
            (m.x[k, i] *
             pyo.exp(m.CapAm - m.CapBm / (m.T[k, i] + m.CapCm)) *
             m.CapBm / (m.T[k, i] + m.CapCm) ** 2 +
             (1 - m.x[k, i]) *
             pyo.exp(m.CapAn - m.CapBn / (m.T[k, i] + m.CapCn)) *
             m.CapBn / (m.T[k, i] + m.CapCn) ** 2
             )

    # Summation equation with tray efficiency(alpha)
    def _gy_rule(m, k, i):
        if k == 1:
            return m.p[1] * m.y[1, i] == m.x[1, i] * pm[1, i]
        elif k == m.Ntray:
            return pyo.Constraint.Skip
        else:
            return m.y[k, i] == \
                m.alpha[k] * m.x[k, i] * pm[k, i] / m.p[k] + \
                (1 - m.alpha[k]) * m.y[k-1, i]

    def _Vm_bound_rule(m,k, i):
        return (0, Vm[k, i], None)

    def _Mv_bound_rule(m,k, i):
        if k == 1:
            return (8.5 + 1e-06, Mv[k, i], 1e+04)
        elif k == m.Ntray:
            return (0.17 + 1e-06, Mv[k, i], 1e+04)
        else:
            return (0.155 + 1e-06, Mv[k, i], 1e+04)

    def _L_bound_rule(m,k, i):
        return (0, L[k, i], 1.0e3)

    m.de_M = pyo.Constraint(m.tray, m.time, rule=_de_M_rule)
    m.de_x = pyo.Constraint(m.tray, m.time, rule=_de_x_rule)

    m.gh = pyo.Constraint(m.tray, m.time, rule=_gh_rule)
    m.dp = pyo.Constraint(m.tray, m.time, rule=_dp_rule)
    m.lTdot = pyo.Constraint(m.tray, m.time, rule=_lTdot_rule)
    m.gy = pyo.Constraint(m.tray, m.time, rule=_gy_rule)

    m.Vm_bound = pyo.Constraint(m.tray, m.time, rule=_Vm_bound_rule)
    m.Mv_bound = pyo.Constraint(m.tray, m.time, rule=_Mv_bound_rule)
    m.L_bound = pyo.Constraint(m.tray, m.time, rule=_L_bound_rule)

    return m


def custom_objective(m, options):
    stage_cost = lambda m, t: -(1e2 * (1 - m.x[1,t]) - m.Qr[t] * 1e1 + m.Qc[t])
    return stage_cost


# ---- Testing the Model ----
if __name__ == '__main__':
    horizon = 10
    m = pyo.ConcreteModel()
    m.time = ContinuousSet(bounds=(0, 10))

    print('Testing Model Equations')
    # try:
    m = variables_initialize(m)

    discretizer = pyo.TransformationFactory('dae.collocation')
    discretizer.apply_to(m, ncp=5, nfe=5, wrt=m.time, scheme='LAGRANGE-RADAU')

    m = equations_write(m)
    #     print('Equation Writing Successful')
    # except Exception as e:
    #     print(f'Equation Writing Failed: {e}')
    # Define and add time-indexed constraints safely using a closure-safe factory

    model_display_flag = True
    if True:  # Toggle to False to disable model display
        with open("model_output.txt", "w") as f:
            m.pprint(ostream=f)

    assert isinstance(m, pyo.ConcreteModel)
