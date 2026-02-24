import pyomo.environ as pyo
from pyomo.dae import DerivativeVar, ContinuousSet
from pyomo.core.expr import sqrt

from infNMPC_options import _import_settings


def variables_initialize(m):
    """
    Initialize variables for a two-column ternary distillation system.

    System: Two columns in series separating three components (A, B, C)
    - Column 1: Separates light A from B+C (produces D1 rich in A, B1 feeds Column 2)
    - Column 2: Separates B from C (produces D2 rich in B, B2 rich in C)

    Components: A (light), B (medium), C (heavy)
    Relative volatilities: alpha_A = 2.0, alpha_B = 1.5 (relative to C)
    """
    # ---- Variable Type Indexing ----
    m.MV_index = pyo.Set(initialize=["VB1", "LT1", "VB2", "LT2"])
    m.CV_index = pyo.Set(initialize=["x1[41,1]", "x2[41,2]", "xC"])
    m.DV_index = pyo.Set(initialize=["F", "qF", "pV", "zF[1]", "zF[2]"])

    m.time_display_name = ["Time (h)"]
    m.MV_display_names = ["VB_1", "LT_1", "VB_2", "LT_2"]
    m.CV_display_names = ["x_{D1,A}", "x_{D2,B}", "x_{B2,C}"]

    # ---- System Parameters ----
    m.NT = pyo.Param(initialize=41)  # Number of trays
    m.NC = pyo.Param(initialize=3)   # Number of components
    m.NF = pyo.Param(initialize=21)  # Feed tray location
    m.tray = pyo.Set(initialize=range(1, 42))  # Tray indices 1..41
    m.vapor_trays = pyo.Set(initialize=range(1, 41))  # Vapor flow defined on trays 1..40
    m.liquid_trays = pyo.Set(initialize=range(2, 42))  # Liquid flow defined on trays 2..41
    m.comp_all = pyo.Set(initialize=[1, 2, 3])  # Component indices (1=A, 2=B, 3=C)
    m.comp = pyo.Set(initialize=[1, 2])

    # ---- Relative Volatilities ----
    m.alpha = pyo.Param({1, 2}, initialize={1: 2.0, 2: 1.5})
    m.m_alpha = pyo.Param({1, 2}, {1, 2}, initialize=lambda m, i, j: m.alpha[i] if i == j else 0)  # Relative volatility matrix (alpha[i] if i==j else 0)

    # ---- Francis Weir Formula Parameters ----
    m.Kuf = pyo.Param(initialize=21.65032)  # Constant above feed
    m.Kbf = pyo.Param(initialize=29.65032)  # Constant below feed
    m.Muw = pyo.Param(initialize=0.25)      # Liquid holdup under weir (kmol)

    # ---- Feed Parameters (Disturbances) ----
    # NOTE: F, qF, zF are mutable and time-indexed. They are defined in
    # equations_write() (called after discretization) so that all discretized
    # time points exist when the Params are constructed and receive their
    # initialize= values. Defining them here would leave post-discretization
    # time points uninitialized (Pyomo does not backfill mutable Params).

    # ---- Consistent steady-state initialization ----
    # clip ≈ M - Muw  (for M > Muw, since sqrt((M-Muw)^2+eps) ≈ M-Muw)
    # L = K * clip^1.5  →  M = Muw + (L/K)^(2/3)
    #
    # Column 1: VB1=4.008, LT1=3.437, F=1.4, qF=1 → D1 = VB1-LT1 = 0.571, B1 = F-D1 = 0.829
    # Below-feed liquid (k=2..21): L1_bf = LT1 + qF*F = 3.437 + 1.4 = 4.837
    #   Weir (Kbf): M1_bf = 0.25 + (4.837/29.65)^(2/3) ≈ 0.548
    # Above-feed liquid (k=22..41): L1_uf = LT1 = 3.437
    #   Weir (Kuf): M1_uf = 0.25 + (3.437/21.65)^(2/3) ≈ 0.543
    #
    # Column 2: VB2=2.40367, LT2=2.13827, B1=0.829 (all-liquid feed from Col-1 bottoms)
    # → D2 = VB2-LT2 = 0.265, B2 = B1-D2 = 0.564
    # Below-feed liquid (k=2..21): L2_bf = LT2 + B1 = 2.138 + 0.829 = 2.967
    #   Weir (Kbf): M2_bf = 0.25 + (2.967/29.65)^(2/3) ≈ 0.465
    # Above-feed liquid (k=22..41): L2_uf = LT2 = 2.138
    #   Weir (Kuf): M2_uf = 0.25 + (2.138/21.65)^(2/3) ≈ 0.464
    #
    # VLE at x1=x2=0.4: num=(0.8,0.6), den=0.4*(2-1)+0.4*(1.5-1)+1=1.6 → y=(0.5,0.375)
    # TC = 0.4*353.3 + 0.4*383.8 + 0.2*411.5 = 377.14
    # xC = 1 - x2[1,1] - x2[1,2] = 1 - 0.4 - 0.4 = 0.2
    _M1_bf = 0.5482  # Col-1 holdup below feed  (weir-consistent with L=4.837)
    _M1_uf = 0.5433  # Col-1 holdup above feed  (weir-consistent with L=3.437)
    _L1_bf = 4.837   # Col-1 liquid below feed  (= LT1 + qF*F, material-balance consistent)
    _L1_uf = 3.437   # Col-1 liquid above feed  (= LT1)
    _M2_bf = 0.4655  # Col-2 holdup below feed  (weir-consistent with L=2.967)
    _M2_uf = 0.4636  # Col-2 holdup above feed  (weir-consistent with L=2.138)
    _L2_bf = 2.967   # Col-2 liquid below feed  (= LT2 + B1, material-balance consistent)
    _L2_uf = 2.138   # Col-2 liquid above feed  (= LT2)

    # Component Fraction, Liquid and Vapor Holdup Variables for Column 1
    m.y1 = pyo.Var(m.vapor_trays, m.comp, m.time, initialize = 0.3)
                   # bounds=(0,1), initialize=lambda m, k, j, t: 0.5 if j == 1 else 0.375)
    m.x1 = pyo.Var(m.tray, m.comp, m.time, initialize=0.4) #, bounds=(0, 1))
    m.M1 = pyo.Var(m.tray, m.time, initialize=0.5)
                   # bounds=(0.01, None), initialize=lambda m, k, t: _M1_bf if k <= 21 else _M1_uf)
    m.V1 = pyo.Var(m.vapor_trays, m.time, initialize=1) # initialize=4.008)
    m.L1 = pyo.Var(m.liquid_trays, m.time, initialize=1)
                   # initialize=lambda m, k, t: _L1_bf if k <= 21 else _L1_uf)
    m.D1 = pyo.Var(m.time, initialize=0.57) # 0.571)
    m.B1 = pyo.Var(m.time, initialize=0.83) # 0.829)

    # Slack Variables for Column 1
    m.x1eps = pyo.Var(m.tray, m.comp, m.time, domain=pyo.NonNegativeReals)
    m.M1eps = pyo.Var(m.tray, m.time, domain=pyo.NonNegativeReals)
    m.V1eps = pyo.Var(m.vapor_trays, m.time, domain=pyo.NonNegativeReals)
    m.L1eps = pyo.Var(m.liquid_trays, m.time, domain=pyo.NonNegativeReals)
    m.D1eps = pyo.Var(m.time, domain=pyo.NonNegativeReals)
    m.B1eps = pyo.Var(m.time, domain=pyo.NonNegativeReals)
    m.TC1eps = pyo.Var(m.tray, m.time, domain=pyo.NonNegativeReals)

    # VLE equation split for Column 1
    m.y1_num = pyo.Var(m.vapor_trays, m.comp, m.time, initialize=0.2)
                       # initialize=lambda m, k, j, t: 0.8 if j == 1 else 0.6)
    m.y1_den = pyo.Var(m.vapor_trays, m.comp, m.time, initialize=0.7) # initialize=1.6)
    m.TC1 = pyo.Var(m.tray, m.time, initialize=380) # 377.14)

    # Component Fraction, Liquid and Vapor Holdup Variables for Column 2
    m.y2 = pyo.Var(m.vapor_trays, m.comp, m.time, initialize=0.5)
                   # bounds=(0, 1), initialize=lambda m, k, j, t: 0.5 if j == 1 else 0.375)
    m.x2 = pyo.Var(m.tray, m.comp, m.time, initialize=0.4) # , bounds=(0, 1))
    m.M2 = pyo.Var(m.tray, m.time, initialize=0.5) 
                   # bounds=(0.01, None), initialize=lambda m, k, t: _M2_bf if k <= 21 else _M2_uf)
    m.V2 = pyo.Var(m.vapor_trays, m.time, initialize=1) # initialize=2.40367)
    m.L2 = pyo.Var(m.liquid_trays, m.time, initialize=1)
                   # initialize=lambda m, k, t: _L2_bf if k <= 21 else _L2_uf)
    m.D2 = pyo.Var(m.time, initialize=0.26) # 0.265)
    m.B2 = pyo.Var(m.time, initialize=0.56) # 0.564)

    # Slack Variables for Column 2
    m.x2eps = pyo.Var(m.tray, m.comp, m.time, domain=pyo.NonNegativeReals)
    m.M2eps = pyo.Var(m.tray, m.time, domain=pyo.NonNegativeReals)
    m.V2eps = pyo.Var(m.vapor_trays, m.time, domain=pyo.NonNegativeReals)
    m.L2eps = pyo.Var(m.liquid_trays, m.time, domain=pyo.NonNegativeReals)
    m.D2eps = pyo.Var(m.time, domain=pyo.NonNegativeReals)
    m.B2eps = pyo.Var(m.time, domain=pyo.NonNegativeReals)
    m.TC2eps = pyo.Var(m.tray, m.time, domain=pyo.NonNegativeReals)
    m.xCeps = pyo.Var(m.time, domain=pyo.NonNegativeReals)

    # VLE equation split for Column 2
    m.y2_num = pyo.Var(m.vapor_trays, m.comp, m.time, initialize = 0.3)
                       # initialize=lambda m, k, j, t: 0.8 if j == 1 else 0.6)
    m.y2_den = pyo.Var(m.vapor_trays, m.comp, m.time, initialize=0.8) # 1.6)
    m.TC2 = pyo.Var(m.tray, m.time, initialize=380) # 377.14)

    # ---- Column 1 Derivative Variables ----
    m.x1dot = DerivativeVar(m.x1, initialize=0.0)
    m.M1dot = DerivativeVar(m.M1, initialize=0.0)

    # ---- Column 2 Derivative Variables ----
    m.x2dot = DerivativeVar(m.x2, initialize=0.0)
    m.M2dot = DerivativeVar(m.M2, initialize=0.0)

    # ---- Column 1 Manipulated Variables ----
    m.VB1 = pyo.Var(m.time, initialize=4.008, bounds=(0, 10))    # Boilup
    m.LT1 = pyo.Var(m.time, initialize=3.43656, bounds=(0, 10))    # Reflux

    # ---- Column 2 Manipulated Variables ----
    m.VB2 = pyo.Var(m.time, initialize=2.40367, bounds=(0, 10))    # Boilup
    m.LT2 = pyo.Var(m.time, initialize=2.13827, bounds=(0, 10))    # Reflux

    # ---- Column 2 C Composition CV ----
    m.xC = pyo.Var(m.time, initialize=0.2)  # xC = 1-x2[1,1]-x2[1,2]; at x2=0.4: 1-0.4-0.4=0.2

    # ---- Setpoints for CVs (product purities) ----
    m.setpoints = pyo.Param(m.CV_index, initialize={
        "x1[41,1]": 0.95,   # A purity in D1
        "x2[41,2]": 0.95,   # B purity in D2
        "xC": 0.95          # C purity in B2
    })

    # ---- Initial Values for CVs ----
    m.initial_values = pyo.Param(m.CV_index, initialize={
        "x1[41,1]": 0.90,
        "x2[41,2]": 0.90,
        "xC": 0.90
    })

    return m


# Make max iters one less than restoration and see most violated constraint.


def equations_write(m):
    """
    Define differential equations and constraints for the two-column system.

    Includes:
    - VLE (vapor-liquid equilibrium) with relative volatilities
    - Material balances (total and component) for each tray
    - Vapor flow relations (constant molar overflow)
    - Liquid flow relations (Francis weir formula)
    """

    # ---- Feed Parameters (Disturbances) ----
    # Defined here (after discretization) so all time points in m.time already
    # exist when the Params are constructed, ensuring initialize= covers them all.
    m.F = pyo.Param(m.time, initialize=1.4, mutable=True)     # Feed flow rate
    m.qF = pyo.Param(m.time, initialize=1.0, mutable=True)    # Feed liquid fraction
    m.zF = pyo.Param(m.comp, m.time, initialize=lambda m, i, t: 0.4 if i == 1 else 0.2, mutable=True)  # Feed composition (A, B)

    # ========================================================================
    # COLUMN 1 EQUATIONS
    # ========================================================================

    # ---- Vapor-Liquid Equilibrium (VLE) for Column 1 ----
    # Using multicomponent ideal VLE with relative volatilities
    # y[i,j] = alpha[j] * x[i,j] / (alpha[1]*x[i,1] + alpha[2]*x[i,2] + 1*(1-x[i,1]-x[i,2]))

    # ---- Soft Constraints for Column 1 ----
    def x1lower_rule(m, k, j, t):
        """Lower bound on liquid composition (with slack)"""
        return m.x1[k, j, t] >= 0 - m.x1eps[k, j, t]
    m.x1_lower = pyo.Constraint(m.tray, m.comp, m.time, rule=x1lower_rule)

    def M1lower_rule(m, k, t):
        """Lower bound on holdup (with slack)"""
        return m.M1[k, t] >= 0.25 - m.M1eps[k, t]
    m.M1_lower = pyo.Constraint(m.tray, m.time, rule=M1lower_rule)

    def M1upper_rule(m, k, t):
        """Upper bound on holdup (with slack)"""
        return m.M1[k, t] <= 0.75 + m.M1eps[k, t]
    m.M1_upper = pyo.Constraint(m.tray, m.time, rule=M1upper_rule)

    def V1lower_rule(m, k, t):
        """Lower bound on vapor flow (with slack)"""
        return m.V1[k, t] >= 0 - m.V1eps[k, t]
    m.V1_lower = pyo.Constraint(m.vapor_trays, m.time, rule=V1lower_rule)

    def L1lower_rule(m, k, t):
        """Lower bound on liquid flow (with slack)"""
        return m.L1[k, t] >= 0 - m.L1eps[k, t]
    m.L1_lower = pyo.Constraint(m.liquid_trays, m.time, rule=L1lower_rule)

    def D1lower_rule(m, t):
        """Lower bound on distillate flow (with slack)"""
        return m.D1[t] >= 0 - m.D1eps[t]
    m.D1_lower = pyo.Constraint(m.time, rule=D1lower_rule)

    def B1lower_rule(m, t):
        """Lower bound on bottoms flow (with slack)"""
        return m.B1[t] >= 0 - m.B1eps[t]
    m.B1_lower = pyo.Constraint(m.time, rule=B1lower_rule)

    # ---- Hard Constraints on TC1 (tray temperature) ----
    def TC1lower_rule(m, k, t):
        """Lower bound on tray temperature (with slack)"""
        return m.TC1[k, t] >= 350
    m.TC1_lower = pyo.Constraint(m.tray, m.time, rule=TC1lower_rule)

    def TC1upper_rule(m, k, t):
        """Upper bound on tray temperature (with slack)"""
        return m.TC1[k, t] <= 450
    m.TC1_upper = pyo.Constraint(m.tray, m.time, rule=TC1upper_rule)

    # ---- VLE Equations for Column 1 ----
    def const1_rule(m, k, j, t):
        """VLE equation numerator"""
        return m.y1_num[k, j, t] == m.x1[k, j, t] * m.m_alpha[j,j]
    m.const1 = pyo.Constraint(m.vapor_trays, m.comp, m.time, rule=const1_rule)

    def const2_rule(m, k, j, t):
        """VLE equation denominator"""
        return m.y1_den[k, j, t] == ((m.x1[k, 1, t] * (m.alpha[1] - 1) + m.x1[k, 2, t]*(m.alpha[2] - 1)) + 1)
    m.const2 = pyo.Constraint(m.vapor_trays, m.comp, m.time, rule=const2_rule)

    def const3_rule(m, k, j, t):
        """VLE equation"""
        return m.y1[k, j, t] * m.y1_den[k, j, t] == m.y1_num[k, j, t]
    m.const3 = pyo.Constraint(m.vapor_trays, m.comp, m.time, rule=const3_rule)

    # ---- Vapor flows assuming constant molar flow for Column 1 ----
    def const4_rule(m, k, t):
        """Vapor flow"""
        if 1 <= k <= m.NF - 1:
            return m.V1[k, t] == m.VB1[t] # Below feed, constant vapor flow equal to boilup
        elif m.NF <= k <= m.NT - 1:
            return m.V1[k, t] == m.VB1[t] + (1 - m.qF[t]) * m.F[t] # Above feed, vapor flow increases by vapor portion of feed
        else:
            return pyo.Constraint.Skip
    m.const4 = pyo.Constraint(m.vapor_trays, m.time, rule=const4_rule)
    
    # ---- Liquid flows using Francis Weir formula for Column 1 ----
    def const6_rule(m, k, t):
        """Liquid flow using Francis Weir formula"""
        if 2 <= k <= m.NF:
            # Below feed
            return m.L1[k, t] == m.Kbf * ( ( m.M1[k,t] - m.Muw + sqrt((m.M1[k, t] - m.Muw)**2 + 10**(-8))) / 2 ) **1.5
        elif m.NF + 1 <= k <= 40:
            # Above feed
            return m.L1[k, t] == m.Kuf * ( ( m.M1[k,t] - m.Muw + sqrt((m.M1[k,t] - m.Muw)**2 + 10**(-8))) / 2 ) **1.5
        elif k == 41:
            # Condenser's liquid flow
            return m.L1[k, t] == m.LT1[t]
        else:
            return pyo.Constraint.Skip
    m.const6 = pyo.Constraint(m.liquid_trays, m.time, rule=const6_rule)

    # ---- Material balances for total holdup and component holdup for Column 1 ----
    def const11_rule(m, k, t):
        """Material balances for total holoup"""
        if k == 1:
            # Reboiler
            return m.M1dot[k, t] == m.L1[2,t] - m.V1[1,t] - m.B1[t]
        elif 2 <= k <= m.NF - 1:
            # Below feed
            return m.M1dot[k, t] == m.L1[k + 1, t] - m.L1[k, t] + m.V1[k - 1, t] - m.V1[k, t]
        elif k == m.NF:
            # Feed tray
            return m.M1dot[k, t] == m.L1[k + 1, t] - m.L1[k, t] + m.V1[k - 1, t] - m.V1[k, t] + m.F[t]
        elif m.NF + 1 <= k <= m.NT - 1:
            # Above feed
            return m.M1dot[k, t] == m.L1[k + 1, t] - m.L1[k, t] + m.V1[k - 1, t] - m.V1[k, t] # Same as below feed
        elif k == m.NT:
            # Total condenser
            return m.M1dot[k, t] == m.V1[k - 1, t] - m.LT1[t] - m.D1[t]
        else:
            return pyo.Constraint.Skip
    m.const11 = pyo.Constraint(m.tray, m.time, rule=const11_rule)

    def const13_rule(m, k, j, t):
        """Material balances for component holdup"""
        if k == 1:
            # Reboiler
            return m.x1dot[k, j, t] * m.M1[k, t] + m.M1dot[k, t] * m.x1[k, j, t] == m.L1[2, t] * m.x1[2, j, t] - m.V1[1, t] * m.y1[1, j, t] - m.B1[t] * m.x1[1, j, t]
        elif 2 <= k <= m.NF - 1:
            # Below feed
            return m.x1dot[k, j, t] * m.M1[k, t] + m.M1dot[k, t] * m.x1[k, j, t] == m.L1[k + 1, t] * m.x1[k + 1, j, t] - m.L1[k, t] * m.x1[k, j, t] + m.V1[k - 1, t] * m.y1[k - 1, j, t] - m.V1[k, t] * m.y1[k, j, t]
        elif k == m.NF:
            # Feed tray
            return m.x1dot[k, j, t] * m.M1[k, t] + m.M1dot[k, t] * m.x1[k, j, t] == m.L1[k + 1, t] * m.x1[k + 1, j, t] - m.L1[k, t] * m.x1[k, j, t] + m.V1[k - 1, t] * m.y1[k - 1, j, t] - m.V1[k, t] * m.y1[k, j, t] + m.F[t] * m.zF[j, t]
        elif m.NF + 1 <= k <= m.NT - 1:
            # Above feed
            return m.x1dot[k, j, t] * m.M1[k, t] + m.M1dot[k, t] * m.x1[k, j, t] == m.L1[k + 1, t] * m.x1[k + 1, j, t] - m.L1[k, t] * m.x1[k, j, t] + m.V1[k - 1, t] * m.y1[k - 1, j, t] - m.V1[k, t] * m.y1[k, j, t]
        elif k == m.NT:
            # Total condenser
            return m.x1dot[k, j, t] * m.M1[k, t] + m.M1dot[k, t] * m.x1[k, j, t] == m.V1[k - 1, t] * m.y1[k - 1, j, t] - m.L1[k, t] * m.x1[k, j, t] - m.D1[t] * m.x1[k, j, t]
        else:
            return pyo.Constraint.Skip
    m.const13 = pyo.Constraint(m.tray, m.comp, m.time, rule=const13_rule)

    # ========================================================================
    # COLUMN 2 EQUATIONS
    # ========================================================================

    # ---- Vapor-Liquid Equilibrium (VLE) for Column 2 ----
    # Using multicomponent ideal VLE with relative volatilities
    # y[i,j] = alpha[j] * x[i,j] / (alpha[1]*x[i,1] + alpha[2]*x[i,2] + 1*(1-x[i,1]-x[i,2]))

    # ---- Soft Constraints for Column 2 ----
    def x2lower_rule(m, k, j, t):
        """Lower bound on liquid composition (with slack)"""
        return m.x2[k, j, t] >= 0 - m.x2eps[k, j, t]
    m.x2_lower = pyo.Constraint(m.tray, m.comp, m.time, rule=x2lower_rule)

    def M2lower_rule(m, k, t):
        """Lower bound on holdup (with slack)"""
        return m.M2[k, t] >= 0.25 - m.M2eps[k, t]
    m.M2_lower = pyo.Constraint(m.tray, m.time, rule=M2lower_rule)

    def M2upper_rule(m, k, t):
        """Upper bound on holdup (with slack)"""
        return m.M2[k, t] <= 0.75 + m.M2eps[k, t]
    m.M2_upper = pyo.Constraint(m.tray, m.time, rule=M2upper_rule)

    def V2lower_rule(m, k, t):
        """Lower bound on vapor flow (with slack)"""
        return m.V2[k, t] >= 0 - m.V2eps[k, t]
    m.V2_lower = pyo.Constraint(m.vapor_trays, m.time, rule=V2lower_rule)

    def L2lower_rule(m, k, t):
        """Lower bound on liquid flow (with slack)"""
        return m.L2[k, t] >= 0 - m.L2eps[k, t]
    m.L2_lower = pyo.Constraint(m.liquid_trays, m.time, rule=L2lower_rule)

    def D2lower_rule(m, t):
        """Lower bound on distillate flow (with slack)"""
        return m.D2[t] >= 0 - m.D2eps[t]
    m.D2_lower = pyo.Constraint(m.time, rule=D2lower_rule)

    def B2lower_rule(m, t):
        """Lower bound on bottoms flow (with slack)"""
        return m.B2[t] >= 0 - m.B2eps[t]
    m.B2_lower = pyo.Constraint(m.time, rule=B2lower_rule)

    # ---- Hard Constraints on TC2 (tray temperature) ----
    def TC2lower_rule(m, k, t):
        """Lower bound on tray temperature (with slack)"""
        return m.TC2[k, t] >= 350
    m.TC2_lower = pyo.Constraint(m.tray, m.time, rule=TC2lower_rule)

    def TC2upper_rule(m, k, t):
        """Upper bound on tray temperature (with slack)"""
        return m.TC2[k, t] <= 450
    m.TC2_upper = pyo.Constraint(m.tray, m.time, rule=TC2upper_rule)

    # ---- VLE Equations for Column 2 ----
    def const1b_rule(m, k, j, t):
        """VLE equation numerator"""
        return m.y2_num[k, j, t] == m.x2[k, j, t] * m.m_alpha[j,j]
    m.const1b = pyo.Constraint(m.vapor_trays, m.comp, m.time, rule=const1b_rule)

    def const2b_rule(m, k, j, t):
        """VLE equation denominator"""
        return m.y2_den[k, j, t] == ((m.x2[k, 1, t] * (m.alpha[1] - 1) + m.x2[k, 2, t]*(m.alpha[2] - 1)) + 1)
    m.const2b = pyo.Constraint(m.vapor_trays, m.comp, m.time, rule=const2b_rule)

    def const3b_rule(m, k, j, t):
        """VLE equation"""
        return m.y2[k, j, t] * m.y2_den[k, j, t] == m.y2_num[k, j, t]
    m.const3b = pyo.Constraint(m.vapor_trays, m.comp, m.time, rule=const3b_rule)

    # ---- Vapor flows assuming constant molar flow for Column 2 ----
    # Col-2 feed (B1) is all liquid (qF2=1), so V2=VB2 everywhere.
    # Consistent with AMPL const25 (k<NF) and const26 (k>=NF) both setting V2=VB2.
    def const4b_rule(m, k, t):
        """Vapor flow = boilup for all trays (all-liquid feed to Col 2)"""
        return m.V2[k, t] == m.VB2[t]
    m.const4b = pyo.Constraint(m.vapor_trays, m.time, rule=const4b_rule)
    
    # ---- Liquid flows using Francis Weir formula for Column 2 ----
    def const6b_rule(m, k, t):
        """Liquid flow using Francis Weir formula"""
        if 2 <= k <= m.NF:
            # Below feed
            return m.L2[k, t] == m.Kbf * ( ( m.M2[k,t] - m.Muw + sqrt((m.M2[k, t] - m.Muw)**2 + 10**(-8))) / 2 ) **1.5
        elif m.NF + 1 <= k <= 40:
            # Above feed
            return m.L2[k, t] == m.Kuf * ( ( m.M2[k,t] - m.Muw + sqrt((m.M2[k,t] - m.Muw)**2 + 10**(-8))) / 2 ) **1.5
        elif k == 41:
            # Condenser's liquid flow
            return m.L2[k, t] == m.LT2[t]
        else:
            return pyo.Constraint.Skip
    m.const6b = pyo.Constraint(m.liquid_trays, m.time, rule=const6b_rule)

    # ---- Material balances for total holdup and component holdup for Column 2 ----
    def const11b_rule(m, k, t):
        """Material balances for total holoup"""
        if k == 1:
            # Reboiler
            return m.M2dot[k, t] == m.L2[2,t] - m.V2[1,t] - m.B2[t]
        elif 2 <= k <= m.NF - 1:
            # Below feed
            return m.M2dot[k, t] == m.L2[k + 1, t] - m.L2[k, t] + m.V2[k - 1, t] - m.V2[k, t]
        elif k == m.NF:
            # Feed tray
            return m.M2dot[k, t] == m.L2[k + 1, t] - m.L2[k, t] + m.V2[k - 1, t] - m.V2[k, t] + m.B1[t]
        elif m.NF + 1 <= k <= m.NT - 1:
            # Above feed
            return m.M2dot[k, t] == m.L2[k + 1, t] - m.L2[k, t] + m.V2[k - 1, t] - m.V2[k, t] # Same as below feed
        elif k == m.NT:
            # Total condenser
            return m.M2dot[k, t] == m.V2[k - 1, t] - m.LT2[t] - m.D2[t]
        else:
            return pyo.Constraint.Skip
    m.const11b = pyo.Constraint(m.tray, m.time, rule=const11b_rule)

    def const13b_rule(m, k, j, t):
        """Material balances for component holdup"""
        if k == 1:
            # Reboiler
            return m.x2dot[k, j, t] * m.M2[k, t] + m.M2dot[k, t] * m.x2[k, j, t] == m.L2[2, t] * m.x2[2, j, t] - m.V2[1, t] * m.y2[1, j, t] - m.B2[t] * m.x2[1, j, t]
        elif 2 <= k <= m.NF - 1:
            # Below feed
            return m.x2dot[k, j, t] * m.M2[k, t] + m.M2dot[k, t] * m.x2[k, j, t] == m.L2[k + 1, t] * m.x2[k + 1, j, t] - m.L2[k, t] * m.x2[k, j, t] + m.V2[k - 1, t] * m.y2[k - 1, j, t] - m.V2[k, t] * m.y2[k, j, t]
        elif k == m.NF:
            # Feed tray
            return m.x2dot[k, j, t] * m.M2[k, t] + m.M2dot[k, t] * m.x2[k, j, t] == m.L2[k + 1, t] * m.x2[k + 1, j, t] - m.L2[k, t] * m.x2[k, j, t] + m.V2[k - 1, t] * m.y2[k - 1, j, t] - m.V2[k, t] * m.y2[k, j, t] + m.B1[t] * m.x1[1,j,t]
        elif m.NF + 1 <= k <= m.NT - 1:
            # Above feed
            return m.x2dot[k, j, t] * m.M2[k, t] + m.M2dot[k ,t]*m.x2 [k,j,t]==m.L2 [k+1,t]*m.x2 [k+1,j,t]-m.L2 [k,t]*m.x2 [k,j,t]+m.V2 [k-1,t]*m.y2 [k-1,j,t]-m.V2 [k,t]*m.y2 [k,j,t]
        elif k == m.NT:
            # Total condenser
            return m.x2dot[k, j, t] * m.M2[k, t] + m.M2dot[k, t] * m.x2[k, j, t] == m.V2[k - 1, t] * m.y2[k - 1, j, t] - m.L2[k, t] * m.x2[k, j, t] - m.D2[t] * m.x2[k, j, t]
        else:
            return pyo.Constraint.Skip
    m.const13b = pyo.Constraint(m.tray, m.comp, m.time, rule=const13b_rule)

    def const43_rule(m, k, t):
        return m.TC1[k, t] == m.x1[k, 1, t] * 353.3 + m.x1[k, 2, t] * 383.8 + (1 - m.x1[k, 1, t] - m.x1[k, 2, t]) * 411.5
    m.const43 = pyo.Constraint(m.tray, m.time, rule=const43_rule)

    def const44_rule(m, k, t):
        return m.TC2[k, t] == m.x2[k, 1, t] * 353.3 + m.x2[k, 2, t] * 383.8 + (1 - m.x2[k, 1, t] - m.x2[k, 2, t]) * 411.5 
    m.const44 = pyo.Constraint(m.tray, m.time, rule=const44_rule)

    def c_bottom_rule(m, t):
        return m.xC[t] == 1 - (m.x2[1, 1, t] + m.x2[1, 2, t])
    m.c_bottom = pyo.Constraint(m.time, rule=c_bottom_rule)

    def const45_rule(m, t):
        return m.x1[41,1,t] >= 0.95 - m.x1eps[41,1,t]
    m.const45 = pyo.Constraint(m.time, rule=const45_rule)
    
    def const46_rule(m, t):
        return m.x2[41,2,t] >= 0.95 - m.x1eps[41,2,t]
    m.const46 = pyo.Constraint(m.time, rule=const46_rule)

    def const47_rule(m, t):
        return m.xC[t] >= 0.95 - m.xCeps[t]
    m.const47 = pyo.Constraint(m.time, rule=const47_rule)

    return m


def custom_objective(m, options):
    """
    Economic objective function for ternary distillation.

    Minimizes: Feed cost + Energy cost - Product revenue
    Cost = pF*F + pV*(VB1 + VB2) - pA*D1 - pB*D2 - pC*B2

    where pA=1, pB=2, pC=1, pF=1, pV is energy price
    """
    # Price parameters
    pF = 1.0   # Feed price
    m.pV = pyo.Param(initialize=1.0)   # Energy price (assumed, can be adjusted)
    pA = 1.0   # Light component A price
    pB = 2.0   # Medium component B price (higher value)
    pC = 1.0   # Heavy component C price
    rho = 1.0e4  # Penalty weight for slack variables (consistent with AMPL rho=1e4)

    def stage_cost(m, t):
        economic_cost = (
            pF * m.F[t] +
            m.pV * (m.VB1[t] + m.VB2[t]) -
            pA * m.D1[t] -
            pB * m.D2[t] -
            pC * m.B2[t]
        )
        penalty = rho * (
            sum(m.x1eps[k, j, t] for k in m.tray for j in m.comp) +
            sum(m.M1eps[k, t] for k in m.tray) +
            sum(m.V1eps[k, t] for k in m.vapor_trays) +
            sum(m.L1eps[k, t] for k in m.liquid_trays) +
            m.D1eps[t] + m.B1eps[t] +
            sum(m.TC1eps[k, t] for k in m.tray) +
            sum(m.x2eps[k, j, t] for k in m.tray for j in m.comp) +
            sum(m.M2eps[k, t] for k in m.tray) +
            sum(m.V2eps[k, t] for k in m.vapor_trays) +
            sum(m.L2eps[k, t] for k in m.liquid_trays) +
            m.D2eps[t] + m.B2eps[t] +
            sum(m.TC2eps[k, t] for k in m.tray) +
            m.xCeps[t]
        )
        return economic_cost + penalty

    return stage_cost


# ---- Testing the Model ----
if __name__ == '__main__':
    m = pyo.ConcreteModel()
    m.time = ContinuousSet(bounds=(0, 10))

    print('Testing Ternary Distillation Model')
    try:
        m = variables_initialize(m)
        print('Variable Initialization Successful')

        discretizer = pyo.TransformationFactory('dae.collocation')
        discretizer.apply_to(m, ncp=3, nfe=2, wrt=m.time, scheme='LAGRANGE-RADAU')
        print('Discretization Successful')

        m = equations_write(m)
        print('Equation Writing Successful')
    except Exception as e:
        print(f'Model Setup Failed: {e}')
        import traceback
        traceback.print_exc()

    model_display_flag = True
    if model_display_flag:
        with open("ternary_model_output.txt", "w") as f:
            m.pprint(ostream=f)

    assert isinstance(m, pyo.ConcreteModel)
    print('Model Test Completed Successfully')

    def test_steady_state():
        from make_model import _make_steady_state_model
        from make_model import _solve_steady_state_model
        options = _import_settings()
        m_ss = pyo.ConcreteModel()
        m_ss = _make_steady_state_model(m_ss, options)

        # ---- Square-system solve switch ----
        # True:  Fix MVs + all slacks to 0, deactivate inequality soft constraints.
        #        IPOPT sees only equality constraints → truly square system.
        #        Converges → model equations are primal feasible.
        #        Fails     → primal infeasibility in the equality equations.
        # False: Full optimization with MVs free and slacks active (normal operation).
        square_system_solve = False

        # Names of all slack variable groups and their corresponding soft constraints
        _slack_vars = [
            'x1eps', 'M1eps', 'V1eps', 'L1eps', 'D1eps', 'B1eps', 'TC1eps',
            'x2eps', 'M2eps', 'V2eps', 'L2eps', 'D2eps', 'B2eps', 'TC2eps',
            'xCeps',
        ]
        _ineq_constraints = [
            'x1_lower', 'M1_lower', 'M1_upper', 'V1_lower', 'L1_lower', 'D1_lower', 'B1_lower',
            'TC1_lower', 'TC1_upper',
            'x2_lower', 'M2_lower', 'M2_upper', 'V2_lower', 'L2_lower', 'D2_lower', 'B2_lower',
            'TC2_lower', 'TC2_upper', 'const45', 'const46', 'const47'
        ]

        if square_system_solve:
            # 1. Fix all MVs to their initialized values
            fixed_mvs = []
            fixed_slacks = []
            for mv_name in m_ss.MV_index:
                try:
                    mv_var = getattr(m_ss, mv_name)
                except AttributeError:
                    print(f"  WARNING: MV '{mv_name}' not found, skipping")
                    continue
                for idx in mv_var:
                    mv_var[idx].fix(pyo.value(mv_var[idx]))
                    fixed_mvs.append(mv_var[idx])

            # 2. Fix all slack variables to 0 (eliminate DOF from penalty terms)
            for sname in _slack_vars:
                sv = getattr(m_ss, sname, None)
                if sv is None:
                    continue
                for idx in sv:
                    sv[idx].fix(0.0)
                    fixed_slacks.append(sv[idx])

            # 3. Deactivate inequality soft constraints (now trivially satisfied)
            deactivated = []
            for cname in _ineq_constraints:
                con = getattr(m_ss, cname, None)
                if con is not None:
                    con.deactivate()
                    deactivated.append(cname)

            print(f"\n=== Square-system solve: fixed {len(fixed_mvs+fixed_slacks)} vars, "
                  f"deactivated {len(deactivated)} inequality constraints ===")

            steady_state_data = _solve_steady_state_model(m_ss, None, options)

            if True:
                with open("base_solution.txt", "w") as f:
                    m_ss.pprint(ostream=f)

            test_flag = True

            # Restore everything
            for v in fixed_mvs:
                # 2: if v.parent_component() is m_ss.VB2 or v.parent_component() is m_ss.LT2:
                # 3: if v.parent_component() is m_ss.VB2 or v.parent_component() is m_ss.LT2 or v.parent_component() is m_ss.VB1:
                # 4: if v.parent_component() is m_ss.VB1:
                # 5: if v.parent_component() is m_ss.VB1 or v.parent_component() is m_ss.LT1:
                # 6: if v.parent_component() is m_ss.LT1:
                # 7: if v.parent_component() is m_ss.VB2 or v.parent_component() is m_ss.LT2 or v.parent_component() is m_ss.LT1:
                if True or test_flag:
                    v.unfix()
                else:
                    continue
            for v in fixed_slacks:

                if test_flag:
                    v.unfix()

                else:

                    comp = v.parent_component()
                    idx = v.index()

                    if (
                        comp is m_ss.M1eps
                        or comp is m_ss.M2eps
                        or comp is m_ss.xCeps
                        or (comp is m_ss.x1eps and idx[0] == 41 and idx[1] == 1)
                        or (comp is m_ss.x2eps and idx[0] == 41 and idx[1] == 2)
                    ):
                        v.unfix()
                    else:
                        continue

            for cname in deactivated:
                if test_flag:
                    getattr(m_ss, cname).activate()
                else:
                    if cname in ['M1_upper', 'M2_upper', 'M1_lower', 'M2_lower', 'const45', 'const46', 'const47']:
                        getattr(m_ss, cname).activate()
                    else:
                        continue
            # for cname in deactivated:
            #     getattr(m_ss, cname).activate()

            steady_state_data = _solve_steady_state_model(m_ss, None, options)

            if True:
                with open("optimized_solution.txt", "w") as f:
                    m_ss.pprint(ostream=f)

        else:
            m_ss = _solve_steady_state_model(m_ss, None, options)

    test_steady_state()
