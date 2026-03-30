"""
eNMPC ternary distillation two-column system.

Two distillation columns in series separating three components (A, B, C).
- Column 1: Separates light A from B+C (produces D1 rich in A, B1 feeds Column 2)
- Column 2: Separates B from C (produces D2 rich in B, B2 rich in C)

Components: A (light), B (medium), C (heavy)
Relative volatilities: alpha_A = 2.0, alpha_B = 1.5 (relative to C)

States: x1[tray,comp] (liquid composition), M1[tray] (liquid holdup) per tray for Col 1;
        x2[tray,comp], M2[tray] for Col 2  (246 state variables total)
MVs:    VB1, LT1, D1, B1 (Column 1); VB2, LT2, D2, B2 (Column 2)
CVs:    xD1A (A purity in D1), xD2B (B purity in D2), xC (C purity in B2)
Objective: minimize feed + energy cost minus product revenue
           cost = pF*F + pV*(VB1+VB2) - pA*D1 - pB*D2 - pC*B2
"""
import pyomo.environ as pyo
from pyomo.dae import DerivativeVar
from pyomo.core.expr import sqrt


def variables_initialize(m):
    """
    Declare all Pyomo variables, sets, and parameters for the two-column system.

    Parameters
    ----------
    m : pyo.ConcreteModel
        m.time is already set by the framework — do NOT re-declare it.

    Returns
    -------
    pyo.ConcreteModel
    """
    print('Initializing Variables')

    # ---- Time Units ----
    m.time_display_name = ["Time (h)"]

    # ---- Variable Type Indexing ----
    m.MV_index = pyo.Set(initialize=["VB1", "LT1", "D1", "B1", "VB2", "LT2", "D2", "B2"])
    m.MV_display_names = ["VB_1", "LT_1", "D_1", "B_1", "VB_2", "LT_2", "D_2", "B_2"]
    m.CV_index = pyo.Set(initialize=["xD1A", "xD2B", "xC"])
    m.CV_display_names = ["x_{D1,A}", "x_{D2,B}", "x_{B2,C}"]
    m.slack_index = pyo.Set(initialize=["xD1A_eps", "xD2B_eps", "xC_eps", "M1_eps", "M2_eps"])

    # ---- System Parameters ----
    m.NT = pyo.Param(initialize=41)   # Number of trays
    m.NC = pyo.Param(initialize=3)    # Number of components
    m.NF = pyo.Param(initialize=21)   # Feed tray location
    m.tray = pyo.Set(initialize=range(1, 42))  # Tray indices 1..41
    m.comp = pyo.Set(initialize=[1, 2])         # Component indices (1=A, 2=B; C is calculated)

    # ---- Relative Volatilities ----
    m.alpha = pyo.Param(m.comp, initialize={1: 2.0, 2: 1.5})

    # ---- Francis Weir Formula Parameters ----
    m.Kuf = pyo.Param(initialize=21.65032)  # Liquid flow constant above feed
    m.Kbf = pyo.Param(initialize=29.65032)  # Liquid flow constant below feed
    m.Muw = pyo.Param(initialize=0.25)      # Liquid holdup under weir (kmol)

    # ---- Feed Parameters (Disturbances) ----
    m.F = pyo.Param(initialize=1.4, mutable=True)                          # Feed flow rate
    m.qF = pyo.Param(initialize=1.0, mutable=True)                        # Feed liquid fraction
    m.zF = pyo.Param(m.comp, initialize={1: 0.4, 2: 0.2}, mutable=True)   # Feed composition

    # ---- Column 1 State Variables ----
    def _x1_init(m, k, j, i):
        if k <= 16:
            return 0.999 * k / 41
        else:
            return 0.36 + (0.98 - 0.36) / (41 - 21 - 1) * (k - 21 - 1)

    def _M1_init(m, k, i):
        if k == 1:
            return 2 * 0.5
        elif 2 <= k <= 21:
            return 0.4 + (0.5 - 0.4) / (21 - 2) * (k - 2)
        else:
            return 0.4 + (0.6 - 0.4) / (41 - 21 - 1) * (k - 21 - 1)

    m.x1 = pyo.Var(m.tray, m.comp, m.time, initialize=_x1_init, bounds=(0, 1))
    m.M1 = pyo.Var(m.tray, m.time, initialize=_M1_init, bounds=(0, 3.0))
    m.x1dot = DerivativeVar(m.x1, initialize=0.0)
    m.M1dot = DerivativeVar(m.M1, initialize=0.0)

    # ---- Column 2 State Variables ----
    def _x2_init(m, k, j, i):
        if k <= 16:
            return 0.5 * k / 41
        else:
            return 0.3 + (0.7 - 0.3) / (41 - 21 - 1) * (k - 21 - 1)

    def _M2_init(m, k, i):
        if k == 1:
            return 2 * 0.5
        elif 2 <= k <= 21:
            return 0.4 + (0.5 - 0.4) / (21 - 2) * (k - 2)
        else:
            return 0.4 + (0.6 - 0.4) / (41 - 21 - 1) * (k - 21 - 1)

    m.x2 = pyo.Var(m.tray, m.comp, m.time, initialize=_x2_init, bounds=(0, 1))
    m.M2 = pyo.Var(m.tray, m.time, initialize=_M2_init, bounds=(0, 3.0))
    m.x2dot = DerivativeVar(m.x2, initialize=0.0)
    m.M2dot = DerivativeVar(m.M2, initialize=0.0)

    # ---- Column 1 Manipulated Variables ----
    m.VB1 = pyo.Var(m.time, initialize=4.008, bounds=(0, 10))   # Boilup
    m.LT1 = pyo.Var(m.time, initialize=3.437, bounds=(0, 10))   # Reflux
    m.D1 = pyo.Var(m.time, initialize=0.57, bounds=(0, 5))      # Distillate
    m.B1 = pyo.Var(m.time, initialize=0.83, bounds=(0, 5))      # Bottoms

    # ---- Column 2 Manipulated Variables ----
    m.VB2 = pyo.Var(m.time, initialize=2.404, bounds=(0, 10))   # Boilup
    m.LT2 = pyo.Var(m.time, initialize=2.138, bounds=(0, 10))   # Reflux
    m.D2 = pyo.Var(m.time, initialize=0.26, bounds=(0, 5))      # Distillate
    m.B2 = pyo.Var(m.time, initialize=0.56, bounds=(0, 5))      # Bottoms

    # ---- Algebraic Variables (Vapor flows — constant molar overflow) ----
    m.V1 = pyo.Var(m.tray, m.time, initialize=4.0, bounds=(0, 20))
    m.V2 = pyo.Var(m.tray, m.time, initialize=2.4, bounds=(0, 20))

    # ---- Controlled Variable Algebraic Variables ----
    # Declared as Var (not Expression) so the framework can load/extract data by key.
    # Equality constraints linking these to the underlying state are added in equations_write.
    m.xD1A = pyo.Var(m.time, initialize=0.9, bounds=(0, 1))  # A purity in Col 1 distillate
    m.xD2B = pyo.Var(m.time, initialize=0.9, bounds=(0, 1))  # B purity in Col 2 distillate
    m.xC = pyo.Var(m.time, initialize=0.9, bounds=(0, 1))    # C purity in Col 2 bottoms

    # ---- Slack Variables (soft constraints, ncp=1 via slack_index) ----
    # CV purity lower bound slacks (one per time point, piecewise-constant)
    m.xD1A_eps = pyo.Var(m.time, initialize=0.0, bounds=(0, None))
    m.xD2B_eps = pyo.Var(m.time, initialize=0.0, bounds=(0, None))
    m.xC_eps   = pyo.Var(m.time, initialize=0.0, bounds=(0, None))
    # Holdup bound slacks
    m.M1_eps = pyo.Var(m.tray, m.time, initialize=0.0, bounds=(0, None))
    m.M2_eps = pyo.Var(m.tray, m.time, initialize=0.0, bounds=(0, None))

    # ---- Setpoints for CVs ----
    m.setpoints = pyo.Param(m.CV_index, initialize={
        "xD1A": 0.95,
        "xD2B": 0.95,
        "xC": 0.95,
    })

    # ---- Initial Values for CVs ----
    m.initial_values = pyo.Param(m.CV_index, initialize={
        "xD1A": 0.90,
        "xD2B": 0.90,
        "xC": 0.90,
    })

    return m


def equations_write(m):
    """
    Add differential equations and algebraic constraints for the two-column system.

    Called AFTER discretization — m.time already contains all collocation points.
    Constraint rules that reference DerivativeVars guard against deleted non-collocation
    entries using ``if (k, t) not in m.M1dot: return pyo.Constraint.Skip``.

    Parameters
    ----------
    m : pyo.ConcreteModel

    Returns
    -------
    pyo.ConcreteModel
    """
    print('Writing Model Equations')

    def clipping(x):
        """Smooth approximation of max(x, 0)."""
        return 0.5 * (x + sqrt(x**2 + 1e-8))

    # ---- CV Equality Constraints ----
    # Link algebraic CV Vars to the underlying tray compositions.
    m.xD1A_def = pyo.Constraint(
        m.time, rule=lambda m, t: m.xD1A[t] == m.x1[41, 1, t]
    )
    m.xD2B_def = pyo.Constraint(
        m.time, rule=lambda m, t: m.xD2B[t] == m.x2[41, 2, t]
    )
    m.xC_def = pyo.Constraint(
        m.time, rule=lambda m, t: m.xC[t] == 1 - m.x2[1, 1, t] - m.x2[1, 2, t]
    )

    # ========================================================================
    # COLUMN 1 EQUATIONS
    # ========================================================================

    # ---- Vapor-Liquid Equilibrium (VLE) for Column 1 ----
    # Multicomponent ideal VLE with relative volatilities (trays 1..40; 41 is total condenser)
    m.y1 = {}
    for i in m.time:
        for k in range(1, 41):
            for j in m.comp:
                m.y1[k, j, i] = (m.alpha[j] * m.x1[k, j, i]) / (
                    m.alpha[1] * m.x1[k, 1, i] +
                    m.alpha[2] * m.x1[k, 2, i] +
                    (1 - m.x1[k, 1, i] - m.x1[k, 2, i])
                )

    # ---- Vapor Flows Column 1 (constant molar overflow) ----
    NF = pyo.value(m.NF)

    def V1_below_feed_rule(m, k, t):
        if 1 <= k <= NF - 1:
            return m.V1[k, t] == m.VB1[t]
        return pyo.Constraint.Skip

    def V1_above_feed_rule(m, k, t):
        if NF <= k <= 40:
            return m.V1[k, t] == m.VB1[t] + (1 - m.qF) * m.F
        return pyo.Constraint.Skip

    m.V1_below_feed = pyo.Constraint(m.tray, m.time, rule=V1_below_feed_rule)
    m.V1_above_feed = pyo.Constraint(m.tray, m.time, rule=V1_above_feed_rule)

    # ---- Liquid Flows Column 1 (Francis weir formula) ----
    m.L1 = {}
    for i in m.time:
        for k in m.tray:
            if 2 <= k <= NF:
                m.L1[k, i] = m.Kbf * clipping(m.M1[k, i] - m.Muw)**1.5
            elif NF + 1 <= k <= 40:
                m.L1[k, i] = m.Kuf * clipping(m.M1[k, i] - m.Muw)**1.5
            elif k == 41:
                m.L1[k, i] = m.LT1[i]

    # ---- Component Holdup Rates Column 1 (algebraic) ----
    m.Mx1dot = {}
    for i in m.time:
        for k in m.tray:
            for j in m.comp:
                if k == 1:
                    # Reboiler
                    m.Mx1dot[k, j, i] = (
                        m.L1[2, i] * m.x1[2, j, i] -
                        m.V1[1, i] * m.y1[1, j, i] -
                        m.B1[i] * m.x1[1, j, i]
                    )
                elif 2 <= k <= NF - 1:
                    # Below feed
                    m.Mx1dot[k, j, i] = (
                        m.L1[k+1, i] * m.x1[k+1, j, i] - m.L1[k, i] * m.x1[k, j, i] +
                        m.V1[k-1, i] * m.y1[k-1, j, i] - m.V1[k, i] * m.y1[k, j, i]
                    )
                elif k == NF:
                    # Feed tray
                    m.Mx1dot[k, j, i] = (
                        m.L1[k+1, i] * m.x1[k+1, j, i] - m.L1[k, i] * m.x1[k, j, i] +
                        m.V1[k-1, i] * m.y1[k-1, j, i] - m.V1[k, i] * m.y1[k, j, i] +
                        m.F * m.zF[j]
                    )
                elif NF + 1 <= k <= 40:
                    # Above feed
                    m.Mx1dot[k, j, i] = (
                        m.L1[k+1, i] * m.x1[k+1, j, i] - m.L1[k, i] * m.x1[k, j, i] +
                        m.V1[k-1, i] * m.y1[k-1, j, i] - m.V1[k, i] * m.y1[k, j, i]
                    )
                elif k == 41:
                    # Total condenser (vapor from tray 40 condenses; no VLE at condenser)
                    m.Mx1dot[k, j, i] = (
                        m.V1[40, i] * m.y1[40, j, i] -
                        m.L1[41, i] * m.x1[41, j, i] -
                        m.D1[i] * m.x1[41, j, i]
                    )

    # ---- Total Material Balance Column 1 ----
    def M1_balance_rule(m, k, t):
        if (k, t) not in m.M1dot:
            return pyo.Constraint.Skip
        if k == 1:
            return m.M1dot[1, t] == m.L1[2, t] - m.V1[1, t] - m.B1[t]
        elif 2 <= k <= NF - 1:
            return m.M1dot[k, t] == m.L1[k+1, t] - m.L1[k, t] + m.V1[k-1, t] - m.V1[k, t]
        elif k == NF:
            return m.M1dot[k, t] == m.L1[k+1, t] - m.L1[k, t] + m.V1[k-1, t] - m.V1[k, t] + m.F
        elif NF + 1 <= k <= 40:
            return m.M1dot[k, t] == m.L1[k+1, t] - m.L1[k, t] + m.V1[k-1, t] - m.V1[k, t]
        elif k == 41:
            return m.M1dot[41, t] == m.V1[40, t] - m.LT1[t] - m.D1[t]
        return pyo.Constraint.Skip

    m.M1_balance = pyo.Constraint(m.tray, m.time, rule=M1_balance_rule)

    # ---- Component Material Balance Column 1 ----
    def x1dot_rule(m, k, j, t):
        if (k, j, t) not in m.x1dot:
            return pyo.Constraint.Skip
        return (
            m.x1dot[k, j, t] * m.M1[k, t] ==
            m.Mx1dot[k, j, t] - m.x1[k, j, t] * m.M1dot[k, t]
        )

    m.x1dot_balance = pyo.Constraint(m.tray, m.comp, m.time, rule=x1dot_rule)

    # ========================================================================
    # COLUMN 2 EQUATIONS
    # ========================================================================

    # ---- VLE for Column 2 (trays 1..40) ----
    m.y2 = {}
    for i in m.time:
        for k in range(1, 41):
            for j in m.comp:
                m.y2[k, j, i] = (m.alpha[j] * m.x2[k, j, i]) / (
                    m.alpha[1] * m.x2[k, 1, i] +
                    m.alpha[2] * m.x2[k, 2, i] +
                    (1 - m.x2[k, 1, i] - m.x2[k, 2, i])
                )

    # ---- Vapor Flows Column 2 (no feed, constant throughout) ----
    def V2_flow_rule(m, k, t):
        if 1 <= k <= 40:
            return m.V2[k, t] == m.VB2[t]
        return pyo.Constraint.Skip

    m.V2_flow = pyo.Constraint(m.tray, m.time, rule=V2_flow_rule)

    # ---- Liquid Flows Column 2 (Francis weir formula) ----
    m.L2 = {}
    for i in m.time:
        for k in m.tray:
            if 2 <= k <= NF:
                m.L2[k, i] = m.Kbf * clipping(m.M2[k, i] - m.Muw)**1.5
            elif NF + 1 <= k <= 40:
                m.L2[k, i] = m.Kuf * clipping(m.M2[k, i] - m.Muw)**1.5
            elif k == 41:
                m.L2[k, i] = m.LT2[i]

    # ---- Component Holdup Rates Column 2 (algebraic) ----
    m.Mx2dot = {}
    for i in m.time:
        for k in m.tray:
            for j in m.comp:
                if k == 1:
                    # Reboiler
                    m.Mx2dot[k, j, i] = (
                        m.L2[2, i] * m.x2[2, j, i] -
                        m.V2[1, i] * m.y2[1, j, i] -
                        m.B2[i] * m.x2[1, j, i]
                    )
                elif 2 <= k <= NF - 1:
                    m.Mx2dot[k, j, i] = (
                        m.L2[k+1, i] * m.x2[k+1, j, i] - m.L2[k, i] * m.x2[k, j, i] +
                        m.V2[k-1, i] * m.y2[k-1, j, i] - m.V2[k, i] * m.y2[k, j, i]
                    )
                elif k == NF:
                    # Feed tray receives B1 (bottoms from Column 1)
                    m.Mx2dot[k, j, i] = (
                        m.L2[k+1, i] * m.x2[k+1, j, i] - m.L2[k, i] * m.x2[k, j, i] +
                        m.V2[k-1, i] * m.y2[k-1, j, i] - m.V2[k, i] * m.y2[k, j, i] +
                        m.B1[i] * m.x1[1, j, i]
                    )
                elif NF + 1 <= k <= 40:
                    m.Mx2dot[k, j, i] = (
                        m.L2[k+1, i] * m.x2[k+1, j, i] - m.L2[k, i] * m.x2[k, j, i] +
                        m.V2[k-1, i] * m.y2[k-1, j, i] - m.V2[k, i] * m.y2[k, j, i]
                    )
                elif k == 41:
                    # Total condenser
                    m.Mx2dot[k, j, i] = (
                        m.V2[40, i] * m.y2[40, j, i] -
                        m.L2[41, i] * m.x2[41, j, i] -
                        m.D2[i] * m.x2[41, j, i]
                    )

    # ---- Total Material Balance Column 2 ----
    def M2_balance_rule(m, k, t):
        if (k, t) not in m.M2dot:
            return pyo.Constraint.Skip
        if k == 1:
            return m.M2dot[1, t] == m.L2[2, t] - m.V2[1, t] - m.B2[t]
        elif 2 <= k <= NF - 1:
            return m.M2dot[k, t] == m.L2[k+1, t] - m.L2[k, t] + m.V2[k-1, t] - m.V2[k, t]
        elif k == NF:
            return m.M2dot[k, t] == m.L2[k+1, t] - m.L2[k, t] + m.V2[k-1, t] - m.V2[k, t] + m.B1[t]
        elif NF + 1 <= k <= 40:
            return m.M2dot[k, t] == m.L2[k+1, t] - m.L2[k, t] + m.V2[k-1, t] - m.V2[k, t]
        elif k == 41:
            return m.M2dot[41, t] == m.V2[40, t] - m.LT2[t] - m.D2[t]
        return pyo.Constraint.Skip

    m.M2_balance = pyo.Constraint(m.tray, m.time, rule=M2_balance_rule)

    # ---- Component Material Balance Column 2 ----
    def x2dot_rule(m, k, j, t):
        if (k, j, t) not in m.x2dot:
            return pyo.Constraint.Skip
        return (
            m.x2dot[k, j, t] * m.M2[k, t] ==
            m.Mx2dot[k, j, t] - m.x2[k, j, t] * m.M2dot[k, t]
        )

    m.x2dot_balance = pyo.Constraint(m.tray, m.comp, m.time, rule=x2dot_rule)

    # ========================================================================
    # SOFT CONSTRAINTS (slack variables)
    # ========================================================================

    # CV purity lower bounds (soft) — mirrors AMPL const45/46/47
    m.xD1A_lb = pyo.Constraint(
        m.time, rule=lambda m, t: m.xD1A[t] >= m.setpoints['xD1A'] - m.xD1A_eps[t]
    )
    m.xD2B_lb = pyo.Constraint(
        m.time, rule=lambda m, t: m.xD2B[t] >= m.setpoints['xD2B'] - m.xD2B_eps[t]
    )
    m.xC_lb = pyo.Constraint(
        m.time, rule=lambda m, t: m.xC[t] >= m.setpoints['xC'] - m.xC_eps[t]
    )

    # Holdup lower bounds (soft) — mirrors AMPL M1lower/M2lower
    def M1_lower_rule(m, k, t):
        lb = 0.5 if k == 1 else 0.25
        return m.M1[k, t] >= lb - m.M1_eps[k, t]

    def M1_upper_rule(m, k, t):
        ub = 2.0 if k == 1 else 0.75
        return m.M1[k, t] <= ub + m.M1_eps[k, t]

    def M2_lower_rule(m, k, t):
        lb = 0.5 if k == 1 else 0.25
        return m.M2[k, t] >= lb - m.M2_eps[k, t]

    def M2_upper_rule(m, k, t):
        ub = 2.0 if k == 1 else 0.75
        return m.M2[k, t] <= ub + m.M2_eps[k, t]

    m.M1_lower = pyo.Constraint(m.tray, m.time, rule=M1_lower_rule)
    m.M1_upper = pyo.Constraint(m.tray, m.time, rule=M1_upper_rule)
    m.M2_lower = pyo.Constraint(m.tray, m.time, rule=M2_lower_rule)
    m.M2_upper = pyo.Constraint(m.tray, m.time, rule=M2_upper_rule)

    return m


def custom_objective(m, options):
    """
    Economic stage cost for ternary distillation with slack penalties.

    Objective = economic cost + L1 penalty on soft constraint violations
        econ    = pF*F + pV*(VB1+VB2) - pA*D1 - pB*D2 - pC*B2
        penalty = rho * (CV purity slacks + holdup bound slacks)
    """
    pF  = 1.0    # Feed price
    pV  = 1.0    # Energy price (per unit of boilup)
    pA  = 1.0    # Light component A product price
    pB  = 2.0    # Medium component B product price (higher value)
    pC  = 1.0    # Heavy component C product price
    rho = 1.0e4  # Slack penalty weight (matches AMPL rho)

    def stage_cost(m, t):
        econ = (
            pF * m.F +
            pV * (m.VB1[t] + m.VB2[t]) -
            pA * m.D1[t] -
            pB * m.D2[t] -
            pC * m.B2[t]
        )
        cv_penalty = rho * (m.xD1A_eps[t] + m.xD2B_eps[t] + m.xC_eps[t])
        M_penalty = rho * (
            sum(m.M1_eps[k, t] for k in m.tray) +
            sum(m.M2_eps[k, t] for k in m.tray)
        )
        return econ + cv_penalty + M_penalty

    return stage_cost


# def default_options():
#     """
#     Model-level default Options overrides.

#     These are applied by Options.for_model_module() before any caller-supplied
#     keyword arguments.
#     """
#     return {
#         'num_horizons': 10,
#         'sampling_time': 1,
#         'nfe_finite': 2,
#         'ncp_finite': 3,
#         'infinite_horizon': True,
#         'nfe_infinite': 3,
#         'ncp_infinite': 3,
#         'tee_flag': True,
#         'endpoint_constraints': True,
#         'custom_objective': True,
#         'input_suppression': True,
#         'input_suppression_factor': 1.0e3,
#         # Weights order: [xD1A, xD2B, xC, VB1, LT1, D1, B1, VB2, LT2, D2, B2]
#         'stage_cost_weights': [1.0e4, 1.0e4, 1.0e4, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
#         'gamma': 0.05,
#         'beta': 1.0,
#         'save_data': True,
#         'save_figure': True,
#     }
