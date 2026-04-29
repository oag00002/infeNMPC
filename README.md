# infeNMPC

[![Tests](https://github.com/oag00002/infeNMPC/actions/workflows/tests.yml/badge.svg)](https://github.com/oag00002/infeNMPC/actions/workflows/tests.yml)

Infinite-horizon economic Nonlinear Model Predictive Control (NMPC) framework built on [Pyomo](https://www.pyomo.org/) and [IDAES](https://idaes.org/).

The distinguishing feature is **infinite-horizon MPC**: the cost-to-go beyond the explicit prediction horizon is approximated by a time-compressed terminal cost integral using the transformation τ = tanh(γt). This gives the controller foresight over the infinite future without solving an infinite-dimensional problem. An optional **Lyapunov stability constraint** can be enforced to provide a certificate of closed-loop stability.

---

## Installation

### Prerequisites

This project requires a forked version of Pyomo (for `clean_model='delete'` DAE support) and IPOPT with the MA57 linear solver (provided by IDAES).

**Step 1 — Install IDAES and IPOPT:**

```bash
pip install idaes-pse tqdm
idaes get-extensions        # downloads IPOPT + MA57 binaries
```

**Step 2 — Install the forked Pyomo:**

```bash
pip install git+https://github.com/oag00002/pyomo.git@mpc-rewrite
```

**Step 3 — Install infeNMPC:**

```bash
git clone https://github.com/oag00002/infeNMPC.git
cd infeNMPC
pip install -e .
```

---

## Quick Start

Write (or choose) a model file that implements `variables_initialize` and `equations_write`, then run a closed-loop simulation:

```python
import my_model                          # your model.py
from infeNMPC import Options, mpc_loop

options = Options.for_model_module(
    my_model,
    num_horizons=100,
    sampling_time=1.0,
    nfe_finite=2,
    ncp_finite=1,
    nfe_infinite=5,
    ncp_infinite=1,
    infinite_horizon=True,
    terminal_constraint_type='none',
    objective='economic',
    stage_cost_weights=[1, 1, 1/600],
    beta=1.0,
    lyap_flag=True,
    lyap_delta=0.01,
    save_data=True,
    save_figure=True,
)
mpc_loop(options)
```

Results are saved to `Results/{model_name}/{date}/{time}/` and include `io_data.csv`, `sim_data.csv`, trajectory plots, and a self-contained `plot_results.py` script.

---

## Writing a Model

Every model file must implement two functions:

```python
def variables_initialize(m):
    """Declare variables, sets, and parameters on the Pyomo model m."""
    m.MV_index = pyo.Set(initialize=["u"])      # manipulated variables
    m.CV_index = pyo.Set(initialize=["x"])      # controlled variables
    m.x  = pyo.Var(m.time, initialize=0.5)
    m.dx = DerivativeVar(m.x, wrt=m.time)
    m.u  = pyo.Var(m.time, bounds=(0, 1), initialize=0.5)
    m.setpoints      = pyo.Param(m.CV_index, initialize={"x": 1.0})
    m.initial_values = pyo.Param(m.CV_index, initialize={"x": 0.0})
    return m

def equations_write(m):
    """Add ODE constraints (called after discretization)."""
    m.ode = pyo.Constraint(m.time, rule=lambda m, t: m.dx[t] == m.u[t] - m.x[t])
    return m

def custom_objective(m, options):
    """Return a stage cost function (required for objective='economic')."""
    return lambda m, t: (m.x[t] - 1.0)**2 + 0.01 * m.u[t]**2
```

**Important:** do not declare `m.time` — the framework injects it before calling `variables_initialize`.

Validate your model with the built-in checker:

```bash
check-model path/to/my_model.py
```

---

## Key Options

| Option | Default | Description |
|--------|---------|-------------|
| `infinite_horizon` | `True` | Use infinite-horizon approximation |
| `nfe_finite` / `ncp_finite` | 2 / 3 | Finite block discretization |
| `nfe_infinite` / `ncp_infinite` | 3 / 3 | Infinite block discretization |
| `objective` | `'economic'` | `'economic'` (custom stage cost) or `'tracking'` (quadratic) |
| `terminal_constraint_type` | `'hard'` | `'hard'`, `'soft'`, or `'none'` |
| `lyap_flag` | `False` | Enable Lyapunov stability constraint |
| `lyap_constraint_type` | `'hard'` | `'hard'`, `'soft'`, or `'none'` |
| `lyap_delta` | `0.01` | Required fractional Lyapunov decrease per step |
| `revive_run` | `0` | Retry attempts on near-convergent IPOPT failures |

See `infeNMPC/infNMPC_options.py` for the full option list with docstrings.

---

## Examples

Working examples are in `examples/`:

```
examples/
  lyap_flag_examples/
    enmpc_cstr/                   <- canonical working example (run this to test)
    enmpc_ternary_distillation/   <- two-column ternary distillation (246 states)
  standard/
    enmpc_cstr/
    binary_distillation/
    nonisothermal_cstr/
    pendulum/
```

Run the canonical CSTR example from the repo root:

```bash
python examples/lyap_flag_examples/enmpc_cstr/run.py
```

---

## Running Tests

```bash
pytest tests/ -q
```

The test suite (121 tests) covers Options configuration, the model validator, collocation tools, solver utilities, Pyomo model structure invariants (t=0 constraint deletion, infinite-block endpoint deletion, Lyapunov infrastructure), and end-to-end MPC loop convergence.
