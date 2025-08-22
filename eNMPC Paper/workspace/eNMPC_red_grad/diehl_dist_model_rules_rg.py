from pyomo.environ import (
    Constraint,
    Expression,
    sqrt,
    exp,
    )
__author__ = 'Kuan-Han Lin'  #: Feb 2022

# mass balances
def rg_de_M_rule(b, i, k):
    m = b.parent_block()
    if i == m.Nfe.at(-1):
        return Constraint.Skip
    else:
        if 1 < k < m.Ntray:
            return (m.V[i, k - 1] - m.V[i, k] + m.L[i, k + 1] - m.L[i, k] + m.feed[k])
        elif k == 1:
            return (m.L[i, 2] - m.L[i, 1] - m.V[i, 1])
        elif k == m.Ntray:
            return (m.V[i, m.Ntray - 1] - m.L[i, m.Ntray] - m.D[i])


def rg_for_dummy_x_rule(b, i, k):
    m = b.parent_block()
    if i == m.Nfe.at(-1):
        return Constraint.Skip
    else:
        if 1 < k < m.Ntray:
            return -(m.dummy_x[i, k] * m.M[i, k] - \
                    (m.V[i, k - 1] * (m.y[i, k - 1] - m.x[i, k]) +
                     m.L[i, k + 1] * (m.x[i, k + 1] - m.x[i, k]) -
                     m.V[i, k] * (m.y[i, k] - m.x[i, k]) +
                     m.feed[k] * (m.xf - m.x[i, k])))
        elif k == 1:
            return -(m.dummy_x[i, 1]*m.M[i, 1] - \
                    (m.L[i, 2] * (m.x[i, 2] - m.x[i, 1]) -
                     m.V[i, 1] * (m.y[i, 1] - m.x[i, 1])))
        elif k == m.Ntray:
            return -(m.dummy_x[i, m.Ntray]*m.M[i, m.Ntray] - \
                    (m.V[i, m.Ntray - 1] * (m.y[i, m.Ntray - 1] - m.x[i, m.Ntray])))


def rg_de_x_rule(b, i, k):
    m = b.parent_block()
    if i == m.Nfe.at(-1):
        return Constraint.Skip
    else:
        if 1 < k < m.Ntray:
            return m.dummy_x[i, k]
        elif k == 1:
            return m.dummy_x[i, 1]
        elif k == m.Ntray:
            return m.dummy_x[i, m.Ntray]


def rg_hrc_rule(b, i):
    m = b.parent_block()
    if i == m.Nfe.at(-1):
        return Constraint.Skip
    else:
        return m.D[i]*m.Rec[i+m.h] - m.L[i, m.Ntray]


# Energy balance
def rg_gh_rule(b, i, k):
    m = b.parent_block()
    if i == m.Nfe.at(-1):
        return Constraint.Skip
    else:
        if 1 < k < m.Ntray:
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
            ) - (m.V[i, k-1] * (m.hv[i, k-1] - m.hl[i, k]) +
                 m.L[i, k+1] * (m.hl[i, k+1] - m.hl[i, k]) -
                 m.V[i, k] * (m.hv[i, k] - m.hl[i, k]) +
                 m.feed[k] * (m.hf - m.hl[i, k]))

        elif k == 1:
            return m.M[i, 1] * (m.dummy_x[i, 1] * ((m.hlm0 - m.hln0) * m.T[i, 1]**3 + (m.hlma - m.hlna)*m.T[i, 1]**2 +
                                                 (m.hlmb - m.hlnb)*m.T[i, 1] +
                                                 m.hlmc - m.hlnc)
                                + m.Tdot[i, 1] * (3 * m.hln0 * m.T[i, 1]**2 + 2 * m.hlna * m.T[i, 1] + m.hlnb +
                                                   m.x[i, 1] *
                                                   (3 * (m.hlm0 - m.hln0) * m.T[i, 1]**2 + 2*(m.hlma - m.hlna) * m.T[i, 1] +
                                                    m.hlmb - m.hlnb)
                                                   )
                                ) - \
                   (m.L[i, 2] * (m.hl[i, 2] - m.hl[i, 1]) - m.V[i, 1] * (m.hv[i, 1] - m.hl[i, 1]) + 1.0E5*m.Qr[i+m.h])

        elif k == m.Ntray:
            return m.M[i, m.Ntray] * (m.dummy_x[i, m.Ntray] * ((m.hlm0 - m.hln0) * m.T[i, m.Ntray]**3 +
                                                             (m.hlma - m.hlna) * m.T[i, m.Ntray]**2 +
                                                             (m.hlmb - m.hlnb) * m.T[i, m.Ntray] +
                                                             m.hlmc - m.hlnc) + m.Tdot[i, m.Ntray] *
                                      (3 * m.hln0 * m.T[i, m.Ntray]**2 + 2* m.hlna * m.T[i, m.Ntray] + m.hlnb +
                                       m.x[i, m.Ntray] * (3 * (m.hlm0 - m.hln0) * m.T[i, m.Ntray]**2 +
                                                          2 * (m.hlma - m.hlna) * m.T[i, m.Ntray] + m.hlmb - m.hlnb))
                                      ) - \
                   (m.V[i, m.Ntray - 1] * (m.hv[i, m.Ntray - 1] - m.hl[i, m.Ntray]) - m.Qc[i])


def rg_hkl_rule(b, i, k):
    m = b.parent_block()
    if i == m.Nfe.at(-1):
        return Constraint.Skip
    else:
        return m.hl[i, k] - (m.x[i, k]*(m.hlm0*m.T[i, k]**3 + m.hlma * m.T[i, k]**2 + m.hlmb * m.T[i, k] + m.hlmc) +\
               (1 - m.x[i, k])*(m.hln0 * m.T[i, k]**3 + m.hlna*m.T[i, k]**2 + m.hlnb * m.T[i, k] + m.hlnc))


def rg_hkv_rule(b, i, k):
    m = b.parent_block()
    if i == m.Nfe.at(-1):
        return Constraint.Skip
    else:
        if k < m.Ntray:
            return m.hv[i, k] - (m.y[i, k] * (m.hlm0 * m.T[i, k]**3 + m.hlma * m.T[i, k]**2 + m.hlmb * m.T[i, k] +
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
                                 )
        else:
            return Constraint.Skip

def rg_lpself_rule(b, i, k):
    m = b.parent_block()
    if i == m.Nfe.at(-1):
        return Constraint.Skip
    else:
        return m.pm[i, k] - (exp(m.CapAm - m.CapBm/(m.T[i, k] + m.CapCm)))


def rg_lpn_rule(b, i, k):
    m = b.parent_block()
    if i == m.Nfe.at(-1):
        return Constraint.Skip
    else:
        return m.pn[i, k] - (exp(m.CapAn - m.CapBn/(m.T[i, k] + m.CapCn)))


def rg_dp_rule(b, i, k):
    m = b.parent_block()
    if i == m.Nfe.at(-1):
        return Constraint.Skip
    else:
        return m.p[k] - (m.pm[i, k] * m.x[i, k] + (1 - m.x[i, k]) * m.pn[i, k])


def rg_lTdot_rule(b, i, k):
    m = b.parent_block()
    if i == m.Nfe.at(-1):
        return Constraint.Skip
    else:
        return -(m.Tdot[i, k] -\
               (-(m.pm[i, k] - m.pn[i, k]) * m.dummy_x[i, k] / \
                (m.x[i, k] *
                 exp(m.CapAm - m.CapBm/(m.T[i, k] + m.CapCm)) * m.CapBm/(m.T[i, k] + m.CapCm)**2 +
                 (1 - m.x[i, k]) *
                 exp(m.CapAn - m.CapBn/(m.T[i, k] + m.CapCn)) * m.CapBn/(m.T[i, k] + m.CapCn)**2)))


def rg_gy0_rule(b, i):
    m = b.parent_block()
    if i == m.Nfe.at(-1):
        return Constraint.Skip
    else:
        return -(m.p[1] * m.y[i, 1] - m.x[i, 1] * m.pm[i, 1])


def rg_gy_rule(b, i, k):
    m = b.parent_block()
    if i == m.Nfe.at(-1):
        return Constraint.Skip
    else:
        if 1 < k < m.Ntray:
            return -(m.y[i, k] - \
                   (m.alpha[k] * m.x[i, k] * m.pm[i, k] / m.p[k] + (1 - m.alpha[k]) * m.y[i, k - 1]))
        else:
            return Constraint.Skip


def rg_dMV_rule(b, i, k):
    m = b.parent_block()
    if i == m.Nfe.at(-1):
        return Constraint.Skip
    else:
        if 1 < k < m.Ntray:
            return m.Mv[i, k] - m.Vm[i, k] * m.M[i, k]
        elif k == 1:
            return m.Mv1[i] - m.Vm[i, 1] * m.M[i, 1]
        elif k == m.Ntray:
            return m.Mvn[i] - m.Vm[i, m.Ntray] * m.M[i, m.Ntray]


def rg_hyd_rule(b, i, k):
    m = b.parent_block()
    if i == m.Nfe.at(-1):
        return Constraint.Skip
    else:
        if 1 < k < m.Ntray:
            return m.L[i, k] * m.Vm[i, k] - 0.166 * (m.Mv[i, k] - 0.155) ** 1.5

        elif k == 1:
            return m.L[i, 1] * m.Vm[i, 1] - 0.166 * (m.Mv1[i] - 8.5) ** 1.5

        elif k == m.Ntray:
            return m.L[i, m.Ntray] * m.Vm[i, m.Ntray] - 0.166 * (m.Mvn[i] - 0.17) ** 1.5


def rg_dvself_rule(b, i, k):
    m = b.parent_block()
    if i == m.Nfe.at(-1):
        return Constraint.Skip
    else:
        return m.Vm[i, k] - (m.x[i, k] * ((1/2288) * 0.2685**(1 + (1 - m.T[i, k]/512.4)**0.2453)) + \
               (1 - m.x[i, k]) * ((1/1235) * 0.27136**(1 + (1 - m.T[i, k]/536.4)**0.24)))


def construct_rg_model_expr(b):
    m = b.parent_block()

    b.rg_de_M = Expression(m.Nfe, m.tray, rule=rg_de_M_rule)
    b.rg_de_x = Expression(m.Nfe, m.tray, rule=rg_de_x_rule)
    b.rg_for_dummy_x = Expression(m.Nfe, m.tray, rule=rg_for_dummy_x_rule)
    b.rg_hrc = Expression(m.Nfe, rule=rg_hrc_rule)
    b.rg_gh = Expression(m.Nfe, m.tray, rule=rg_gh_rule)
    b.rg_hkl = Expression(m.Nfe, m.tray, rule=rg_hkl_rule)
    b.rg_hkv = Expression(m.Nfe, m.tray, rule=rg_hkv_rule)
    b.rg_lpself = Expression(m.Nfe, m.tray, rule=rg_lpself_rule)
    b.rg_lpn = Expression(m.Nfe, m.tray, rule=rg_lpn_rule)
    b.rg_dp = Expression(m.Nfe, m.tray, rule=rg_dp_rule)
    b.rg_lTdot = Expression(m.Nfe, m.tray, rule=rg_lTdot_rule)
    b.rg_gy0 = Expression(m.Nfe, rule=rg_gy0_rule)
    b.rg_gy = Expression(m.Nfe, m.tray, rule=rg_gy_rule)
    b.rg_dMV = Expression(m.Nfe, m.tray, rule=rg_dMV_rule)
    b.rg_hyd = Expression(m.Nfe, m.tray, rule=rg_hyd_rule)
    b.rg_dvself = Expression(m.Nfe, m.tray, rule=rg_dvself_rule)