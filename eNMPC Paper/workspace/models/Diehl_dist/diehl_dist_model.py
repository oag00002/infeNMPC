from pyomo.environ import (
    ConcreteModel,
    Set,
    Var,
    Param,
    Constraint,
    sqrt,
    exp,
    )
from pyomo.dae import (
    DerivativeVar,
    ContinuousSet,
    )

__author__ = 'David Thierry, Kuan-Han Lin'  #: May 2018
# modified by KHL on Oct 2020

# mass balances
def de_M_rule(m, i, k):
    if i > 0 and 1 < k < m.Ntray:
        return m.Mdot[i, k] == \
               (m.V[i, k - 1] - m.V[i, k] + m.L[i, k + 1] - m.L[i, k] + m.feed[k])
    elif i > 0 and k == 1:
        return m.Mdot[i, 1] == \
               (m.L[i, 2] - m.L[i, 1] - m.V[i, 1])
    elif i > 0 and k == m.Ntray:
        return m.Mdot[i, m.Ntray] == \
               (m.V[i, m.Ntray - 1] - m.L[i, m.Ntray] - m.D[i])
    else:
        return Constraint.Skip


def for_dummy_x_rule(m, i, k):
    if i > 0 and 1 < k < m.Ntray:
        return m.dummy_x[i, k] * m.M[i, k] == \
                (m.V[i, k - 1] * (m.y[i, k - 1] - m.x[i, k]) +
                 m.L[i, k + 1] * (m.x[i, k + 1] - m.x[i, k]) -
                 m.V[i, k] * (m.y[i, k] - m.x[i, k]) +
                 m.feed[k] * (m.xf - m.x[i, k]))
    elif i > 0 and k == 1:
        return m.dummy_x[i, 1]*m.M[i, 1] == \
                (m.L[i, 2] * (m.x[i, 2] - m.x[i, 1]) -
                 m.V[i, 1] * (m.y[i, 1] - m.x[i, 1]))
    elif i > 0 and k == m.Ntray:
        return m.dummy_x[i, m.Ntray]*m.M[i, m.Ntray] == \
                (m.V[i, m.Ntray - 1] * (m.y[i, m.Ntray - 1] - m.x[i, m.Ntray]))
    else:
        return Constraint.Skip


def de_x_rule(m, i, k):
    if i > 0 and 1 < k < m.Ntray:
        return m.xdot[i, k] == m.dummy_x[i, k]
    elif i > 0 and k == 1:
        return m.xdot[i, 1] == m.dummy_x[i, 1]
    elif i > 0 and k == m.Ntray:
        return m.xdot[i, m.Ntray] == m.dummy_x[i, m.Ntray]
    else:
        return Constraint.Skip


def hrc_rule(m, i):
    if i > 0:
        return m.D[i]*m.Rec[i] - m.L[i, m.Ntray] == 0
    else:
        return Constraint.Skip


# Energy balance
def gh_rule(m, i, k):
    if i > 0 and 1 < k < m.Ntray:
        return m.M[i, k] * (
         m.dummy_x[i, k] * (
                (m.hlm0 - m.hln0) * (m.T[i, k]**3) +
                (m.hlma - m.hlna) * (m.T[i, k]**2) +
                (m.hlmb - m.hlnb) * m.T[i, k] + m.hlmc - m.hlnc) +
             m.Tdot[i, k] * (
                3*m.hln0*(m.T[i, k]**2) +
                2*m.hlna * m.T[i, k] + m.hlnb +
                m.x[i, k] *
                (3*(m.hlm0 - m.hln0) * (m.T[i, k]**2) + 2 * (m.hlma - m.hlna) * m.T[i, k] + m.hlmb - m.hlnb))
        ) == (m.V[i, k-1] * (m.hv[i, k-1] - m.hl[i, k]) +
              m.L[i, k+1] * (m.hl[i, k+1] - m.hl[i, k]) -
              m.V[i, k] * (m.hv[i, k] - m.hl[i, k]) +
              m.feed[k] * (m.hf - m.hl[i, k]))

    elif i > 0 and k == 1:
        return m.M[i, 1] * (m.dummy_x[i, 1] * ((m.hlm0 - m.hln0) * m.T[i, 1]**3 + (m.hlma - m.hlna)*m.T[i, 1]**2 +
                                             (m.hlmb - m.hlnb)*m.T[i, 1] +
                                             m.hlmc - m.hlnc)
                            + m.Tdot[i, 1] * (3 * m.hln0 * m.T[i, 1]**2 + 2 * m.hlna * m.T[i, 1] + m.hlnb +
                                               m.x[i, 1] *
                                               (3 * (m.hlm0 - m.hln0) * m.T[i, 1]**2 + 2*(m.hlma - m.hlna) * m.T[i, 1] +
                                                m.hlmb - m.hlnb)
                                               )
                            ) == \
               (m.L[i, 2] * (m.hl[i, 2] - m.hl[i, 1]) - m.V[i, 1] * (m.hv[i, 1] - m.hl[i, 1]) + 1.0E5*m.Qr[i])

    elif i > 0 and k == m.Ntray:
        return m.M[i, m.Ntray] * (m.dummy_x[i, m.Ntray] * ((m.hlm0 - m.hln0) * m.T[i, m.Ntray]**3 +
                                                         (m.hlma - m.hlna) * m.T[i, m.Ntray]**2 +
                                                         (m.hlmb - m.hlnb) * m.T[i, m.Ntray] +
                                                         m.hlmc - m.hlnc) + m.Tdot[i, m.Ntray] *
                                  (3 * m.hln0 * m.T[i, m.Ntray]**2 + 2* m.hlna * m.T[i, m.Ntray] + m.hlnb +
                                   m.x[i, m.Ntray] * (3 * (m.hlm0 - m.hln0) * m.T[i, m.Ntray]**2 +
                                                      2 * (m.hlma - m.hlna) * m.T[i, m.Ntray] + m.hlmb - m.hlnb))
                                  ) == \
               (m.V[i, m.Ntray - 1] * (m.hv[i, m.Ntray - 1] - m.hl[i, m.Ntray]) - m.Qc[i])
    else:
        return Constraint.Skip


def hkl_rule(m, i, k):
    if i > 0:
        return m.hl[i, k] == m.x[i, k]*(m.hlm0*m.T[i, k]**3 + m.hlma * m.T[i, k]**2 + m.hlmb * m.T[i, k] + m.hlmc) +\
               (1 - m.x[i, k])*(m.hln0 * m.T[i, k]**3 + m.hlna*m.T[i, k]**2 + m.hlnb * m.T[i, k] + m.hlnc)
    else:
        return Constraint.Skip


def hkv_rule(m, i, k):
    if i > 0 and k < m.Ntray:
        return m.hv[i, k] == m.y[i, k] * (m.hlm0 * m.T[i, k]**3 + m.hlma * m.T[i, k]**2 + m.hlmb * m.T[i, k] +
                                          m.hlmc +
                                          m.r * m.Tkm * sqrt(1 - (m.p[k]/m.Pkm) * (m.Tkm/m.T[i, k])**3) *
                                          (m.a - m.b * m.T[i, k]/m.Tkm + m.c1 * (m.T[i, k]/m.Tkm)**7 + m.gm *
                                           (m.d - m.l * m.T[i, k]/m.Tkm + m.f*(m.T[i, k]/m.Tkm)**7))
                                          ) + (1 - m.y[i, k]) * (m.hln0 * m.T[i, k]**3 + m.hlna * m.T[i, k]**2 +
                                                                 m.hlnb * m.T[i, k] +
                                                                 m.hlnc +
                                                                 m.r * m.Tkn *
                                                                 sqrt(1 - (m.p[k]/m.Pkn) * (m.Tkn/m.T[i, k])**3) *
                                                                 (m.a - m.b * m.T[i, k]/m.Tkn +
                                                                  m.c1 * (m.T[i, k]/m.Tkn)**7 +
                                                                  m.gn*(m.d - m.l * m.T[i, k]/m.Tkn +
                                                                        m.f* (m.T[i, k]/m.Tkn)**7)
                                                                  )
                                                                 )
    else:
        return Constraint.Skip


def lpself_rule(m, i, k):
    if i > 0:
        return m.pm[i, k] == exp(m.CapAm - m.CapBm/(m.T[i, k] + m.CapCm))
    else:
        return Constraint.Skip


def lpn_rule(m, i, k):
    if i > 0:
        return m.pn[i, k] == exp(m.CapAn - m.CapBn/(m.T[i, k] + m.CapCn))

    else:
        return Constraint.Skip


def dp_rule(m, i, k):
    if i > 0:
        return m.p[k] == m.pm[i, k] * m.x[i, k] + (1 - m.x[i, k]) * m.pn[i, k]
    else:
        return Constraint.Skip


def lTdot_rule(m, i, k):
    if i > 0:
        return m.Tdot[i, k]  ==\
               -(m.pm[i, k] - m.pn[i, k]) * m.dummy_x[i, k] / \
               (m.x[i, k] *
                exp(m.CapAm - m.CapBm/(m.T[i, k] + m.CapCm)) * m.CapBm/(m.T[i, k] + m.CapCm)**2 +
                (1 - m.x[i, k]) *
                exp(m.CapAn - m.CapBn/(m.T[i, k] + m.CapCn)) * m.CapBn/(m.T[i, k] + m.CapCn)**2)
    else:
        return Constraint.Skip


def gy0_rule(m, i):
    if i > 0:
        return m.p[1] * m.y[i, 1] == m.x[i, 1] * m.pm[i, 1]
    else:
        return Constraint.Skip


def gy_rule(m, i, k):
    if i > 0 and 1 < k < m.Ntray:
        return m.y[i, k] == \
               m.alpha[k] * m.x[i, k] * m.pm[i, k] / m.p[k] + (1 - m.alpha[k]) * m.y[i, k - 1]
    else:
        return Constraint.Skip


def dMV_rule(m, i, k):
    if i > 0 and 1 < k < m.Ntray:
        return m.Mv[i, k] == m.Vm[i, k] * m.M[i, k]
    elif i > 0 and k == 1:
        return m.Mv1[i] == m.Vm[i, 1] * m.M[i, 1]
    elif i > 0 and k == m.Ntray:
        return m.Mvn[i] == m.Vm[i, m.Ntray] * m.M[i, m.Ntray]
    else:
        return Constraint.Skip


def hyd_rule(m, i, k):
    if i > 0 and 1 < k < m.Ntray:
        return m.L[i, k] * m.Vm[i, k] == 0.166 * (m.Mv[i, k] - 0.155) ** 1.5

    elif i > 0 and k == 1:
        return m.L[i, 1] * m.Vm[i, 1] == 0.166 * (m.Mv1[i] - 8.5) ** 1.5

    elif i > 0 and k == m.Ntray:
        return m.L[i, m.Ntray] * m.Vm[i, m.Ntray] == 0.166 * (m.Mvn[i] - 0.17) ** 1.5

    else:
        return Constraint.Skip


def dvself_rule(m, i, k):
    if i > 0:
        return m.Vm[i, k] == m.x[i, k] * ((1/2288) * 0.2685**(1 + (1 - m.T[i, k]/512.4)**0.2453)) + \
               (1 - m.x[i, k]) * ((1/1235) * 0.27136**(1 + (1 - m.T[i, k]/536.4)**0.24))
    else:
        return Constraint.Skip

# ---------------------------------------------------------------------------------------------------------------------
def dist_daemodel(nfe = 20, th = 60):
    mod = ConcreteModel()

    mod.t = ContinuousSet(bounds=(0, th*nfe))
    mod.h = Param(initialize = th)
    mod.Ntray = Ntray = 42

    mod.tray = Set(initialize=[i for i in range(1, mod.Ntray + 1)])
    mod.feed = Param(mod.tray,
                      initialize=lambda m, k: 57.5294 if k == 21 else 0.0,
                      mutable=True)

    mod.xf = Param(initialize=0.32, mutable=True)  # feed mole fraction
    mod.hf = Param(initialize=9081.3)  # feed enthalpy

    mod.hlm0 = Param(initialize=2.6786e-04)
    mod.hlma = Param(initialize=-0.14779)
    mod.hlmb = Param(initialize=97.4289)
    mod.hlmc = Param(initialize=-2.1045e04)

    mod.hln0 = Param(initialize=4.0449e-04)
    mod.hlna = Param(initialize=-0.1435)
    mod.hlnb = Param(initialize=121.7981)
    mod.hlnc = Param(initialize=-3.0718e04)

    mod.r = Param(initialize=8.3147)
    mod.a = Param(initialize=6.09648)
    mod.b = Param(initialize=1.28862)
    mod.c1 = Param(initialize=1.016)
    mod.d = Param(initialize=15.6875)
    mod.l = Param(initialize=13.4721)
    mod.f = Param(initialize=2.615)

    mod.gm = Param(initialize=0.557)
    mod.Tkm = Param(initialize=512.6)
    mod.Pkm = Param(initialize=8.096e06)

    mod.gn = Param(initialize=0.612)
    mod.Tkn = Param(initialize=536.7)
    mod.Pkn = Param(initialize=5.166e06)

    mod.CapAm = Param(initialize=23.48)
    mod.CapBm = Param(initialize=3626.6)
    mod.CapCm = Param(initialize=-34.29)

    mod.CapAn = Param(initialize=22.437)
    mod.CapBn = Param(initialize=3166.64)
    mod.CapCn = Param(initialize=-80.15)

    mod.pstrip = Param(initialize=250)
    mod.prect = Param(initialize=190)


    def _p_init(m, k):
        ptray = 9.39e+04
        if k <= 20:
            return _p_init(m, 21) + m.pstrip * (21 - k)
        elif 20 < k < m.Ntray:
            return ptray + m.prect * (m.Ntray - k)
        elif k == m.Ntray:
            return 9.39e+04
    mod.p = Param(mod.tray, initialize=_p_init)

    def _alpha_init(m, i):
        if i <= 21:
            return 0.62
        else:
            return 0.35
    mod.alpha = Param(mod.tray,
                      initialize=lambda m, k: 0.62 if k <= 21 else 0.35)

    # --------------------------------------------------------------------------------------------------------------
    #: First define differential state variables (state: x, ic-Param: x_ic, derivative-Var:xdot
    #: States (differential) section

    def __m_init(m, i, k):
        if k < m.Ntray:
            return 4000.
        elif k == 1:
            return 104340.
        elif k == m.Ntray:
            return 5000.

    #: Liquid hold-up
    mod.M = Var(mod.t, mod.tray, initialize=__m_init)
    #: Mole-fraction
    mod.x = Var(mod.t, mod.tray, initialize=lambda m, i, k: 0.999 * k / m.Ntray)

    #:  Derivative-var
    mod.Mdot = DerivativeVar(mod.M, initialize=0.0)
    mod.xdot = DerivativeVar(mod.x, initialize=0.0)

    mod.dummy_x = Var(mod.t, mod.tray, initialize=0.0)
    # --------------------------------------------------------------------------------------------------------------
    # States (algebraic) section
    # Tray temperature
    mod.T = Var(mod.t, mod.tray,
                initialize=lambda m, i, k: ((370.781 - 335.753) / m.Ntray) * k + 370.781)
    mod.Tdot = Var(mod.t, mod.tray, initialize=1e-05)  #: Not really a der_var

    # saturation pressures
    mod.pm = Var(mod.t, mod.tray, initialize=1e4)
    mod.pn = Var(mod.t, mod.tray, initialize=1e4)

    # Vapor mole flowrate
    mod.V = Var(mod.t, mod.tray, initialize=44.0)

    def _l_init(m, i, k):
        if 2 <= k <= 21:
            return 83.
        elif 22 <= k <= 42:
            return 23
        elif k == 1:
            return 40

    # Liquid mole flowrate
    mod.L = Var(mod.t, mod.tray, initialize=_l_init)

    # Vapor mole frac & diff var
    mod.y = Var(mod.t, mod.tray,
                 initialize=lambda m, i, k: ((0.99 - 0.005) / m.Ntray) * k + 0.005)

    # Liquid enthalpy # enthalpy
    mod.hl = Var(mod.t, mod.tray, initialize=10000.)

    # Liquid enthalpy # enthalpy
    mod.hv = Var(mod.t, mod.tray, initialize=5e+04)
    # Re-boiler & condenser heat
    mod.Qc = Var(mod.t, initialize=1.6e06)
    mod.D = Var(mod.t, initialize=18.33)
    # vol holdups
    mod.Vm = Var(mod.t, mod.tray, initialize=6e-05)

    mod.Mv = Var(mod.t, mod.tray,
                  initialize=lambda m, i, k: 0.23 if 1 < k < m.Ntray else 0.0)#, bounds=(0.1550001, None))
    mod.Mv1 = Var(mod.t, initialize=8.57)#, bounds=(8.50001,None))
    mod.Mvn = Var(mod.t, initialize=0.203)#, bounds=(0.170001,None))

    # --------------------------------------------------------------------------------------------------------------
    #: Controls
    mod.Rec = Var(mod.t, initialize=7.72700925775773761472464684629813E-01)
    mod.Qr = Var(mod.t, initialize=1.78604740940007800236344337463379E+01)

    # --------------------------------------------------------------------------------------------------------------
    #: Constraints for the differential states
    #: Then the ode-Con:de_x, collocation-Con:dvar_t_x, noisy-Expr: noisy_x, cp-Constraint: cp_x, initial-Con: x_icc
    #: Differential equations
    mod.de_M = Constraint(mod.t, mod.tray, rule=de_M_rule)
    mod.de_x = Constraint(mod.t, mod.tray, rule=de_x_rule)
    mod.for_dummy_x = Constraint(mod.t, mod.tray, rule=for_dummy_x_rule)
    # --------------------------------------------------------------------------------------------------------------
    #: Constraint section (algebraic equations)

    mod.hrc = Constraint(mod.t, rule=hrc_rule)
    mod.gh = Constraint(mod.t, mod.tray, rule=gh_rule)
    mod.hkl = Constraint(mod.t, mod.tray, rule=hkl_rule)
    mod.hkv = Constraint(mod.t, mod.tray, rule=hkv_rule)
    mod.lpself = Constraint(mod.t, mod.tray, rule=lpself_rule)
    mod.lpn = Constraint(mod.t, mod.tray, rule=lpn_rule)
    mod.dp = Constraint(mod.t, mod.tray, rule=dp_rule)
    mod.lTdot = Constraint(mod.t, mod.tray, rule=lTdot_rule)
    mod.gy0 = Constraint(mod.t, rule=gy0_rule)
    mod.gy = Constraint(mod.t, mod.tray, rule=gy_rule)
    mod.dMV = Constraint(mod.t, mod.tray, rule=dMV_rule)
    mod.hyd = Constraint(mod.t, mod.tray, rule=hyd_rule)
    mod.dvself = Constraint(mod.t, mod.tray, rule=dvself_rule)

    # if delete_unused_var:
    #     del mod.V[:, 42]
    #     del mod.y[:, 42]
    #     del mod.hv[:, 42]
    #     del mod.Mv[:, 1]
    #     del mod.Mv[:, 42]

    return mod

def set_ics(dist_mod, ic_states):
    for i in mod.tray:
        mod.x_ic[i] = ic_states["x"][i]
        mod.M_ic[i] = ic_states["M"][i]

def initial_w_ics(dist_mod, ic_states):
    for i in mod.t:
        for j in mod.tray:
            mod.x[(i,j)] = ic_states["x"][j]
            mod.M[(i,j)] = ic_states["M"][j]

def create_model_hard_bounds(mod):
    state_bounds = {
                    "M": (1.0, 1e+07),
                    "T": (200, 500),
                    "pm": (1.0, 5e+07),
                    "pn": (1.0, 5e+07),
                    "L": (0.0, 1e+03),
                    "V": (0.0, 1e+03),
                    "x": (0.0, 1.0),
                    "y": (0.0, 1.0),
                    "hl": (1.0, 1e+07),
                    "hv": (1.0, 1e+07),
                    "Qc": (0.0, 1e+08),
                    "D": (0.0, 1e+04),
                    "Vm": (0.0, 1e+04),
                    "Mv": (0.155 + 1e-06, 1e+04),
                    "Mv1": (8.5 + 1e-06, 1e+04),
                    "Mvn": (0.17 + 1e-06, 1e+04)
                    }
    u_bounds = {"Rec": (0000.1, 99.999), "Qr": (0, None)}

    for key, val in {**state_bounds, **u_bounds}.items():
        var = mod.find_component(key)
        var.setlb(val[0])
        var.setub(val[1])
