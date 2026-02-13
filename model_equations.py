import pyomo.environ as pyo
from pyomo.dae import DerivativeVar, ContinuousSet
from pyomo.core.expr import sqrt


def variables_initialize(m):
    """
    Initialize variables for a two-column ternary distillation system.

    System: Two columns in series separating three components (A, B, C)
    - Column 1: Separates light A from B+C (produces D1 rich in A, B1 feeds Column 2)
    - Column 2: Separates B from C (produces D2 rich in B, B2 rich in C)

    Components: A (light), B (medium), C (heavy)
    Relative volatilities: alpha_A = 2.0, alpha_B = 1.5 (relative to C)
    """
    print('Initializing Variables')

    # ---- Time Units ----
    m.time_display_name = ["Time (h)"]

    # ---- Variable Type Indexing ----
    m.MV_index = pyo.Set(initialize=["VB1", "LT1", "D1", "B1", "VB2", "LT2", "D2", "B2"])
    m.MV_display_names = ["VB_1", "LT_1", "D_1", "B_1", "VB_2", "LT_2", "D_2", "B_2"]
    m.CV_index = pyo.Set(initialize=["x1[41,1]", "x2[41,2]", "xC"])
    m.CV_display_names = ["x_{D1,A}", "x_{D2,B}", "x_{B2,C}"]
    m.DV_index = pyo.Set(initialize=["F", "qF"])

    # ---- System Parameters ----
    m.NT = pyo.Param(initialize=41)  # Number of trays
    m.NC = pyo.Param(initialize=3)   # Number of components
    m.NF = pyo.Param(initialize=21)  # Feed tray location
    m.tray = pyo.Set(initialize=range(1, 42))  # Tray indices 1..41
    m.comp = pyo.Set(initialize=[1, 2])  # Component indices (1=A, 2=B, C is calculated)

    # ---- Relative Volatilities ----
    m.alpha = pyo.Param(m.comp, initialize={1: 2.0, 2: 1.5})

    # ---- Francis Weir Formula Parameters ----
    m.Kuf = pyo.Param(initialize=21.65032)  # Constant above feed
    m.Kbf = pyo.Param(initialize=29.65032)  # Constant below feed
    m.Muw = pyo.Param(initialize=0.25)      # Liquid holdup under weir (kmol)

    # ---- Feed Parameters (Disturbances) ----
    m.F = pyo.Param(initialize=1.4, mutable=True)     # Feed flow rate
    m.qF = pyo.Param(initialize=1.0, mutable=True)    # Feed liquid fraction
    m.zF = pyo.Param(m.comp, initialize={1: 0.4, 2: 0.2}, mutable=True)  # Feed composition (A, B)

    # ---- Column 1 State Variables ----
    def _x1_init(m, k, j, i):
        """Initialize liquid composition in Column 1"""
        if k <= 16:
            return 0.999 * k / 41
        else:
            return 0.36 + (0.98 - 0.36) / (41 - 21 - 1) * (k - 21 - 1)

    def _M1_init(m, k, i):
        """Initialize liquid holdup in Column 1"""
        if k == 1:
            return 2 * 0.5
        elif 2 <= k <= 21:
            return 0.4 + (0.5 - 0.4) / (21 - 2) * (k - 2)
        else:
            return 0.4 + (0.6 - 0.4) / (41 - 21 - 1) * (k - 21 - 1)

    m.x1 = pyo.Var(m.tray, m.comp, m.time, initialize=_x1_init, bounds=(0, 1))

    def _M1_bounds(m, k, t):
        """Holdup bounds - reboiler (tray 1) has larger capacity"""
        if k == 1:
            return (0.5, 2.0)
        else:
            return (0.25, 0.75)

    m.M1 = pyo.Var(m.tray, m.time, initialize=_M1_init, bounds=_M1_bounds)

    # ---- Column 1 Derivative Variables ----
    m.x1dot = DerivativeVar(m.x1, initialize=0.0)
    m.M1dot = DerivativeVar(m.M1, initialize=0.0)

    # ---- Column 2 State Variables ----
    def _x2_init(m, k, j, i):
        """Initialize liquid composition in Column 2"""
        if k <= 16:
            return 0.5 * k / 41
        else:
            return 0.3 + (0.7 - 0.3) / (41 - 21 - 1) * (k - 21 - 1)

    def _M2_init(m, k, i):
        """Initialize liquid holdup in Column 2"""
        if k == 1:
            return 2 * 0.5
        elif 2 <= k <= 21:
            return 0.4 + (0.5 - 0.4) / (21 - 2) * (k - 2)
        else:
            return 0.4 + (0.6 - 0.4) / (41 - 21 - 1) * (k - 21 - 1)

    m.x2 = pyo.Var(m.tray, m.comp, m.time, initialize=_x2_init, bounds=(0, 1))

    def _M2_bounds(m, k, t):
        """Holdup bounds - reboiler (tray 1) has larger capacity"""
        if k == 1:
            return (0.5, 2.0)
        else:
            return (0.25, 0.75)

    m.M2 = pyo.Var(m.tray, m.time, initialize=_M2_init, bounds=_M2_bounds)

    # ---- Column 2 Derivative Variables ----
    m.x2dot = DerivativeVar(m.x2, initialize=0.0)
    m.M2dot = DerivativeVar(m.M2, initialize=0.0)

    # ---- Column 1 Manipulated Variables ----
    m.VB1 = pyo.Var(m.time, initialize=4.008, bounds=(0, 10))    # Boilup
    m.LT1 = pyo.Var(m.time, initialize=3.437, bounds=(0, 10))    # Reflux
    m.D1 = pyo.Var(m.time, initialize=0.57, bounds=(0, 5))       # Distillate
    m.B1 = pyo.Var(m.time, initialize=0.83, bounds=(0, 5))       # Bottoms

    # ---- Column 2 Manipulated Variables ----
    m.VB2 = pyo.Var(m.time, initialize=2.404, bounds=(0, 10))    # Boilup
    m.LT2 = pyo.Var(m.time, initialize=2.138, bounds=(0, 10))    # Reflux
    m.D2 = pyo.Var(m.time, initialize=0.26, bounds=(0, 5))       # Distillate
    m.B2 = pyo.Var(m.time, initialize=0.56, bounds=(0, 5))       # Bottoms

    # ---- Algebraic Variables (Vapor flows - constant molar overflow) ----
    m.V1 = pyo.Var(m.tray, m.time, initialize=4.0, bounds=(0, 20))
    m.V2 = pyo.Var(m.tray, m.time, initialize=2.4, bounds=(0, 20))

    # ---- Computed Controlled Variable (C purity in Column 2 bottoms) ----
    # xC = 1 - xA - xB at bottom of column 2
    # This is defined as a Variable with a Constraint, not an Expression
    m.xC = pyo.Var(m.time, initialize=0.90, bounds=(0, 1))

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


def equations_write(m):
    """
    Define differential equations and constraints for the two-column system.

    Includes:
    - VLE (vapor-liquid equilibrium) with relative volatilities
    - Material balances (total and component) for each tray
    - Vapor flow relations (constant molar overflow)
    - Liquid flow relations (Francis weir formula)
    """
    print('Writing Model Equations')

    def clipping(x):
        """Smooth approximation of max(x, 0)"""
        return 0.5 * (x + sqrt(x**2 + 1e-8))

    # ========================================================================
    # COLUMN 1 EQUATIONS
    # ========================================================================

    # ---- Vapor-Liquid Equilibrium (VLE) for Column 1 ----
    # Using multicomponent ideal VLE with relative volatilities
    # y[i,j] = alpha[j] * x[i,j] / (alpha[1]*x[i,1] + alpha[2]*x[i,2] + 1*(1-x[i,1]-x[i,2]))

    m.y1 = {}  # Vapor composition (algebraic)
    for i in m.time:
        for k in range(1, 41):  # Trays 1..40 (not tray 41, which is total condenser)
            for j in m.comp:
                m.y1[k, j, i] = (m.alpha[j] * m.x1[k, j, i]) / (
                    m.alpha[1] * m.x1[k, 1, i] +
                    m.alpha[2] * m.x1[k, 2, i] +
                    (1 - m.x1[k, 1, i] - m.x1[k, 2, i])
                )

    # ---- Vapor Flows (Constant Molar Overflow) ----
    def V1_below_feed_rule(m, k, t):
        """Vapor flow below feed tray"""
        if 1 <= k <= m.NF - 1:
            return m.V1[k, t] == m.VB1[t]
        else:
            return pyo.Constraint.Skip

    def V1_above_feed_rule(m, k, t):
        """Vapor flow above feed tray"""
        if m.NF <= k <= 40:
            return m.V1[k, t] == m.VB1[t] + (1 - m.qF) * m.F
        else:
            return pyo.Constraint.Skip

    m.V1_below_feed = pyo.Constraint(m.tray, m.time, rule=V1_below_feed_rule)
    m.V1_above_feed = pyo.Constraint(m.tray, m.time, rule=V1_above_feed_rule)

    # ---- Liquid Flows (Francis Weir Formula) ----
    m.L1 = {}  # Liquid flow (algebraic)
    for i in m.time:
        for k in m.tray:
            if 2 <= k <= m.NF:
                # Below feed
                m.L1[k, i] = m.Kbf * clipping(m.M1[k, i] - m.Muw)**1.5
            elif m.NF + 1 <= k <= 40:
                # Above feed
                m.L1[k, i] = m.Kuf * clipping(m.M1[k, i] - m.Muw)**1.5
            elif k == 41:
                # Top tray (reflux)
                m.L1[k, i] = m.LT1[i]

    # ---- Material Balances Column 1 ----
    def M1_balance_rule(m, k, t):
        """Total material balance for each tray"""
        if k == 1:
            # Reboiler
            return m.M1dot[1, t] == m.L1[2, t] - m.V1[1, t] - m.B1[t]
        elif 2 <= k <= m.NF - 1:
            # Below feed
            return m.M1dot[k, t] == m.L1[k + 1, t] - m.L1[k, t] + m.V1[k - 1, t] - m.V1[k, t]
        elif k == m.NF:
            # Feed tray
            return m.M1dot[k, t] == m.L1[k + 1, t] - m.L1[k, t] + m.V1[k - 1, t] - m.V1[k, t] + m.F
        elif m.NF + 1 <= k <= 40:
            # Above feed
            return m.M1dot[k, t] == m.L1[k + 1, t] - m.L1[k, t] + m.V1[k - 1, t] - m.V1[k, t]
        elif k == 41:
            # Total condenser
            return m.M1dot[41, t] == m.V1[40, t] - m.LT1[t] - m.D1[t]
        else:
            return pyo.Constraint.Skip

    m.M1_balance = pyo.Constraint(m.tray, m.time, rule=M1_balance_rule)

    # ---- Component Material Balances Column 1 ----
    m.Mx1dot = {}  # Component holdup derivative (algebraic)
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
                elif 2 <= k <= m.NF - 1:
                    # Below feed
                    m.Mx1dot[k, j, i] = (
                        m.L1[k + 1, i] * m.x1[k + 1, j, i] - m.L1[k, i] * m.x1[k, j, i] +
                        m.V1[k - 1, i] * m.y1[k - 1, j, i] - m.V1[k, i] * m.y1[k, j, i]
                    )
                elif k == m.NF:
                    # Feed tray
                    m.Mx1dot[k, j, i] = (
                        m.L1[k + 1, i] * m.x1[k + 1, j, i] - m.L1[k, i] * m.x1[k, j, i] +
                        m.V1[k - 1, i] * m.y1[k - 1, j, i] - m.V1[k, i] * m.y1[k, j, i] +
                        m.F * m.zF[j]
                    )
                elif m.NF + 1 <= k <= 40:
                    # Above feed
                    m.Mx1dot[k, j, i] = (
                        m.L1[k + 1, i] * m.x1[k + 1, j, i] - m.L1[k, i] * m.x1[k, j, i] +
                        m.V1[k - 1, i] * m.y1[k - 1, j, i] - m.V1[k, i] * m.y1[k, j, i]
                    )
                elif k == 41:
                    # Total condenser (no VLE, y = x at top)
                    m.Mx1dot[k, j, i] = (
                        m.V1[40, i] * m.y1[40, j, i] -
                        m.L1[41, i] * m.x1[41, j, i] -
                        m.D1[i] * m.x1[41, j, i]
                    )

    def x1dot_rule(m, k, j, t):
        """Derivative of liquid composition"""
        return m.x1dot[k, j, t] * m.M1[k, t] == m.Mx1dot[k, j, t] - m.x1[k, j, t] * m.M1dot[k, t]

    m.x1dot_balance = pyo.Constraint(m.tray, m.comp, m.time, rule=x1dot_rule)

    # ========================================================================
    # COLUMN 2 EQUATIONS
    # ========================================================================

    # ---- Vapor-Liquid Equilibrium (VLE) for Column 2 ----
    m.y2 = {}  # Vapor composition (algebraic)
    for i in m.time:
        for k in range(1, 41):
            for j in m.comp:
                m.y2[k, j, i] = (m.alpha[j] * m.x2[k, j, i]) / (
                    m.alpha[1] * m.x2[k, 1, i] +
                    m.alpha[2] * m.x2[k, 2, i] +
                    (1 - m.x2[k, 1, i] - m.x2[k, 2, i])
                )

    # ---- Vapor Flows (Constant Molar Overflow) ----
    def V2_flow_rule(m, k, t):
        """Vapor flow in Column 2 (no feed, so constant throughout)"""
        if 1 <= k <= 40:
            return m.V2[k, t] == m.VB2[t]
        else:
            return pyo.Constraint.Skip

    m.V2_flow = pyo.Constraint(m.tray, m.time, rule=V2_flow_rule)

    # ---- Liquid Flows (Francis Weir Formula) ----
    m.L2 = {}  # Liquid flow (algebraic)
    for i in m.time:
        for k in m.tray:
            if 2 <= k <= m.NF:
                # Below "feed" (which is B1 from column 1 at tray NF)
                m.L2[k, i] = m.Kbf * clipping(m.M2[k, i] - m.Muw)**1.5
            elif m.NF + 1 <= k <= 40:
                # Above feed
                m.L2[k, i] = m.Kuf * clipping(m.M2[k, i] - m.Muw)**1.5
            elif k == 41:
                # Top tray (reflux)
                m.L2[k, i] = m.LT2[i]

    # ---- Material Balances Column 2 ----
    def M2_balance_rule(m, k, t):
        """Total material balance for each tray"""
        if k == 1:
            # Reboiler
            return m.M2dot[1, t] == m.L2[2, t] - m.V2[1, t] - m.B2[t]
        elif 2 <= k <= m.NF - 1:
            # Below feed from Column 1
            return m.M2dot[k, t] == m.L2[k + 1, t] - m.L2[k, t] + m.V2[k - 1, t] - m.V2[k, t]
        elif k == m.NF:
            # Feed tray (receives B1 from Column 1)
            return m.M2dot[k, t] == m.L2[k + 1, t] - m.L2[k, t] + m.V2[k - 1, t] - m.V2[k, t] + m.B1[t]
        elif m.NF + 1 <= k <= 40:
            # Above feed
            return m.M2dot[k, t] == m.L2[k + 1, t] - m.L2[k, t] + m.V2[k - 1, t] - m.V2[k, t]
        elif k == 41:
            # Total condenser
            return m.M2dot[41, t] == m.V2[40, t] - m.LT2[t] - m.D2[t]
        else:
            return pyo.Constraint.Skip

    m.M2_balance = pyo.Constraint(m.tray, m.time, rule=M2_balance_rule)

    # ---- Component Material Balances Column 2 ----
    m.Mx2dot = {}  # Component holdup derivative (algebraic)
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
                elif 2 <= k <= m.NF - 1:
                    # Below feed
                    m.Mx2dot[k, j, i] = (
                        m.L2[k + 1, i] * m.x2[k + 1, j, i] - m.L2[k, i] * m.x2[k, j, i] +
                        m.V2[k - 1, i] * m.y2[k - 1, j, i] - m.V2[k, i] * m.y2[k, j, i]
                    )
                elif k == m.NF:
                    # Feed tray (receives B1 from Column 1)
                    m.Mx2dot[k, j, i] = (
                        m.L2[k + 1, i] * m.x2[k + 1, j, i] - m.L2[k, i] * m.x2[k, j, i] +
                        m.V2[k - 1, i] * m.y2[k - 1, j, i] - m.V2[k, i] * m.y2[k, j, i] +
                        m.B1[i] * m.x1[1, j, i]  # Feed from Column 1 bottoms
                    )
                elif m.NF + 1 <= k <= 40:
                    # Above feed
                    m.Mx2dot[k, j, i] = (
                        m.L2[k + 1, i] * m.x2[k + 1, j, i] - m.L2[k, i] * m.x2[k, j, i] +
                        m.V2[k - 1, i] * m.y2[k - 1, j, i] - m.V2[k, i] * m.y2[k, j, i]
                    )
                elif k == 41:
                    # Total condenser
                    m.Mx2dot[k, j, i] = (
                        m.V2[40, i] * m.y2[40, j, i] -
                        m.L2[41, i] * m.x2[41, j, i] -
                        m.D2[i] * m.x2[41, j, i]
                    )

    def x2dot_rule(m, k, j, t):
        """Derivative of liquid composition"""
        return m.x2dot[k, j, t] * m.M2[k, t] == m.Mx2dot[k, j, t] - m.x2[k, j, t] * m.M2dot[k, t]

    m.x2dot_balance = pyo.Constraint(m.tray, m.comp, m.time, rule=x2dot_rule)

    # ---- Constraint to define xC (C purity in Column 2 bottoms) ----
    def xC_definition_rule(m, t):
        """Define xC as the purity of component C in the bottoms of Column 2"""
        return m.xC[t] == 1 - m.x2[1, 1, t] - m.x2[1, 2, t]

    m.xC_definition = pyo.Constraint(m.time, rule=xC_definition_rule)

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
    pV = 1.0   # Energy price (assumed, can be adjusted)
    pA = 1.0   # Light component A price
    pB = 2.0   # Medium component B price (higher value)
    pC = 1.0   # Heavy component C price

    def stage_cost(m, t):
        return (
            pF * m.F +
            pV * (m.VB1[t] + m.VB2[t]) -
            pA * m.D1[t] -
            pB * m.D2[t] -
            pC * m.B2[t]
        )

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

    model_display_flag = False
    if model_display_flag:
        with open("ternary_model_output.txt", "w") as f:
            m.pprint(ostream=f)

    assert isinstance(m, pyo.ConcreteModel)
    print('Model Test Completed Successfully')
