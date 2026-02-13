# Ternary Distillation - Double Column System

## System Description

Two distillation columns in series separating three components (A, B, C):

- **Column 1**: Separates light component A from B+C
  - Top product (D1): Rich in A
  - Bottoms (B1): Feeds Column 2

- **Column 2**: Separates B from C
  - Top product (D2): Rich in B
  - Bottoms (B2): Rich in C

## System Parameters

- **Trays**: 41 per column (NT = 41)
- **Feed location**: Tray 21 (NF = 21)
- **Components**:
  - A (light): α_A = 2.0
  - B (medium): α_B = 1.5
  - C (heavy): α_C = 1.0 (reference)

## States (246 total)

**Column 1** (123 states):
- `x1[tray, comp]`: Liquid composition (41 × 2 = 82 states)
- `M1[tray]`: Liquid holdup (41 states)

**Column 2** (123 states):
- `x2[tray, comp]`: Liquid composition (41 × 2 = 82 states)
- `M2[tray]`: Liquid holdup (41 states)

*Note: Third component C is calculated as: xC = 1 - xA - xB*

## Manipulated Variables (8)

**Column 1**:
- `VB1`: Reboiler boilup rate
- `LT1`: Top reflux rate
- `D1`: Distillate flowrate
- `B1`: Bottoms flowrate

**Column 2**:
- `VB2`: Reboiler boilup rate
- `LT2`: Top reflux rate
- `D2`: Distillate flowrate
- `B2`: Bottoms flowrate

## Controlled Variables (3)

Product purities (implemented as Pyomo Expressions):
- `xD1A`: A purity in Column 1 distillate = `x1[41,1]` (setpoint: 0.95)
- `xD2B`: B purity in Column 2 distillate = `x2[41,2]` (setpoint: 0.95)
- `xC`: C purity in Column 2 bottoms = `1-x2[1,1]-x2[1,2]` (setpoint: 0.95)

## Disturbances (2)

- `F`: Feed flowrate (kmol/h)
- `qF`: Feed liquid fraction (0-1)
- `zF`: Feed composition [zA, zB] (zC = 1 - zA - zB)

## Model Equations

### VLE (Vapor-Liquid Equilibrium)
Ideal multicomponent VLE with relative volatilities:

```
y[i,j] = (α[j] * x[i,j]) / (α[A]*x[i,A] + α[B]*x[i,B] + x[i,C])
```

### Vapor Flow (Constant Molar Overflow)
- **Column 1 below feed**: V1 = VB1
- **Column 1 above feed**: V1 = VB1 + (1-qF)*F
- **Column 2**: V2 = VB2 (constant, no external feed)

### Liquid Flow (Francis Weir Formula)
```
L[i] = K * (M[i] - Muw)^1.5
```

where:
- K = 29.65 (below feed), 21.65 (above feed)
- Muw = 0.25 kmol (holdup under weir)

### Material Balances
For each tray:
```
dM[i]/dt = L[i+1] - L[i] + V[i-1] - V[i] + F_in[i]
d(M*x)[i,j]/dt = L[i+1]*x[i+1,j] - L[i]*x[i,j] + V[i-1]*y[i-1,j] - V[i]*y[i,j] + F_in*z[j]
```

## Economic Objective (custom_objective)

Minimize operating cost:
```
Cost = pF*F + pV*(VB1 + VB2) - pA*D1 - pB*D2 - pC*B2
```

where:
- pF = 1.0 (feed price)
- pV = 1.0 (energy price)
- pA = 1.0 (light product price)
- pB = 2.0 (medium product price - most valuable)
- pC = 1.0 (heavy product price)

## Usage

1. Copy `model_equations.py` and `infNMPC_options.py` to the main directory
2. Run the MPC simulation:
   ```bash
   python run_MPC.py
   ```

## Notes

- Total of 246 differential states makes this a large-scale problem
- Material balance constraints couple the MVs (D1, B1, etc.) through dynamics
- The two columns are coupled: B1 from Column 1 feeds Column 2 at tray NF
- Economic objective favors producing B (pB=2.0) while minimizing energy costs

## References

Based on AMPL model:
- File: `double_col/column2_dynamic_soft.mod`
- Three-component (A, B, C) distillation
- Two columns in series configuration
