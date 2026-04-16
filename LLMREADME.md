# infeNMPC — LLM Orientation Guide

This document is a seed for future AI agents. It describes the architecture, information flow, and implementation details of the **infeNMPC** codebase so that an agent can operate confidently without re-reading every file.

---

## What This Project Is

**infeNMPC** is a Python framework for running closed-loop Nonlinear Model Predictive Control (NMPC) simulations. Its distinguishing feature is support for **infinite-horizon** NMPC, where the cost-to-go beyond the explicit prediction horizon is approximated by a time-transformed terminal cost integral. It is built on **Pyomo** (algebraic modeling) and **IDAES** (process simulation). The solver is always **IPOPT** with MA57.

The two public entry points are:

```python
from infeNMPC import Options, mpc_loop
options = Options.for_model_module(model_module, **kwargs)
mpc_loop(options)
```

---

## Directory Layout

```
infeNMPC/
├── infeNMPC/                        # Main Python package
│   ├── __init__.py                  # Re-exports: Options, Plant, InfiniteHorizonController,
│   │                                #   FiniteHorizonController, mpc_loop
│   ├── __main__.py                  # CLI: python -m infeNMPC --model path/to/model.py
│   ├── infNMPC_options.py           # Options dataclass — all configuration lives here
│   ├── run_MPC.py                   # mpc_loop() — the closed-loop simulation orchestrator
│   ├── controllers.py               # InfiniteHorizonController, FiniteHorizonController
│   ├── plant.py                     # Plant — single-step dynamic simulator
│   ├── make_model.py                # Pyomo model construction (most complex file)
│   ├── model_equations.py           # Model loader: _get_model(), _load_model()
│   ├── check_model.py               # Model validator: static AST + import + dynamic checks
│   ├── data_save_and_plot.py        # CSV export, figure generation, plot_results.py script
│   ├── live_plot.py                 # Live monitoring subprocess (used with safe_run=True)
│   └── tools/
│       ├── collocation_tools.py     # _compute_gamma_from_collocation()
│       ├── debug_tools.py           # _report_constraint_violations()
│       ├── indexing_tools.py        # Variable key helpers, expression builders
│       └── initialization_tools.py  # Progressive sampling-time warm-start
├── examples/
│   ├── standard/                    # Examples without Lyapunov constraint
│   │   ├── enmpc_cstr/              # Economic CSTR: A→B (canonical example, WORKS)
│   │   ├── pendulum/                # Inverted pendulum on cart
│   │   ├── binary_distillation/     # 42-tray distillation column (finite horizon only)
│   │   ├── enmpc_binary_distillation/
│   │   ├── nonisothermal_cstr/
│   │   └── ternary_distillation_double_column/  # Two-column ternary distillation (246 states)
│   ├── lyap_flag_examples/          # Same systems with lyap_flag=True (or optional)
│   │   ├── enmpc_cstr/              # CANONICAL WORKING EXAMPLE — run this to test
│   │   ├── enmpc_binary_distillation/
│   │   └── enmpc_ternary_distillation/  # Two-column ternary distillation with Lyap support
│   └── live_plot_examples/
│       └── enmpc_cstr.py            # Single-file live-monitoring example
└── pyproject.toml                   # Dependencies: pyomo, idaes-pse, tqdm
```

**The only confirmed working test example:** `examples/lyap_flag_examples/enmpc_cstr/run.py`

---

## The Options Dataclass (`infNMPC_options.py`)

`Options` is a `@dataclass` with defaults for every field. Never instantiate it raw — use the classmethods:

```python
Options.for_model_module(module, **overrides)  # preferred: pass live module object
Options.for_model(model_name, **overrides)     # loads from infeNMPC/models/ by name
options.copy(**overrides)                      # immutable-style copy with changes
```

### Key Fields

| Field | Default | Meaning |
|---|---|---|
| `model_module` | `None` | Live Python module (takes priority over `model_name`) |
| `model_name` | `'unknown'` | Dotted import name or display name |
| `num_horizons` | `100` | Number of MPC iterations |
| `sampling_time` | `1.0` | Control update interval (in model time units) |
| `nfe_finite` | `2` | Finite elements in finite-horizon block |
| `ncp_finite` | `3` | Collocation points per finite element (finite block) |
| `infinite_horizon` | `True` | Use infinite-horizon vs. finite-horizon controller |
| `nfe_infinite` | `3` | Finite elements in infinite-horizon block |
| `ncp_infinite` | `3` | Collocation points per element (infinite block) |
| `custom_objective` | `True` | Use model's `custom_objective()` economic cost |
| `terminal_constraint_type` | `'hard'` | `'hard'` = equality at terminal time; `'soft'` = quadratic penalty; `'none'` = disabled |
| `terminal_constraint_variables` | `'cvmv'` | Which variables to constrain: `'cv'`, `'mv'`, or `'cvmv'` |
| `terminal_soft_weight` | `1.0` | Scalar multiplier on soft-constraint penalty (per-var weights from `stage_cost_weights`) |
| `terminal_cost_riemann` | `False` | Riemann-sum approximation of terminal cost |
| `stage_cost_weights` | `[1,1,1/600]` | Per-variable weights for quadratic tracking cost |
| `gamma` | `None` | Time-compression param; **auto-computed if None** |
| `beta` | `1.0` | Weight on terminal cost vs. stage cost |
| `lyap_flag` | `False` | Build Lyapunov infrastructure (`phi_track`, `V_prev`, constraint); constraint form set by `lyap_constraint_type` |
| `lyap_delta` | `0.01` | Required fractional decrease per step (0=inactive, 1=tight) |
| `lyap_constraint_type` | `'hard'` | `'hard'` = strict inequality constraint; `'soft'` = L1-relaxed via `lyap_slack` Var in objective; `'none'` = infrastructure built, no constraint added |
| `lyap_soft_weight` | `1e4` | Weight on L1 slack penalty when `lyap_constraint_type='soft'` |
| `input_suppression` | `False` | Penalize MV increments |
| `input_suppression_factor` | `1.0` | Weight on move-suppression penalty |
| `initialize_with_initial_data` | `False` | Spread IC across all time points |
| `initialization_assist` | `False` | Progressive sampling-time warm-start |
| `initialization_assist_sampling_time_start` | `10.0` | Starting large sampling time |
| `safe_run` | `False` | Write `io_data.csv` after every iteration |
| `save_data` | `True` | Save CSVs at end of run |
| `save_figure` | `True` | Save matplotlib figures |
| `tee_flag` | `False` | Print IPOPT output |
| `disturb_flag` | `False` | Apply random plant disturbances |
| `disturb_distribution` | `'normal'` | Disturbance distribution type |
| `disturb_seeded` | `True` | Use fixed random seed |
| `debug_flag` | `False` | Print most-violated constraints after every controller solve and on failures |
| `dynamic_initial_conditions` | `False` | Use dynamic IC path (fix state vars to `initialize=` values, solve for algebraic vars only) instead of SS-NLP |

---

## The User Model Contract (`model_equations.py`)

Every example has a `model.py` that implements this interface. The framework loads it via `_get_model(options)` which returns `options.model_module` if set, or imports by `options.model_name`.

### Required Functions

```python
def variables_initialize(m):
    """
    Called on a fresh ConcreteModel (or Block). Must add:
      - m.time        (ContinuousSet, already set by framework — DO NOT re-declare)
      - m.MV_index    (pyo.Set of manipulated variable names as strings)
      - m.CV_index    (pyo.Set of controlled variable names as strings)
      - All state variables as pyo.Var(m.time, ...)
      - All MV variables as pyo.Var(m.time, bounds=(...))
      - DerivativeVar for each state variable
      - m.setpoints   (pyo.Param indexed by CV_index, used if not custom_objective)
      - m.initial_values (pyo.Param indexed by CV_index or a subset)
    Optional display metadata:
      - m.time_display_name  = ["Time (units)"]
      - m.CV_display_names   = ["latex_name_1", ...]
      - m.MV_display_names   = ["latex_name_1", ...]
    Returns m.
    """

def equations_write(m):
    """
    Add all differential equations and algebraic constraints.
    The framework calls this AFTER discretization, so m.time already has
    all collocation points populated.

    Do NOT guard against t=0 — the framework deletes all constraint
    entries at t=0 immediately after this call (RADAU: t=0 is not a
    collocation point).  Write rules over m.time without skip logic.
    Returns m.
    """
```

### Optional Functions

```python
def custom_objective(m, options):
    """
    Return a callable: stage_cost_fn(m, t) -> scalar Pyomo expression.
    Used when options.custom_objective=True.
    The steady-state value of this function is subtracted inside the framework
    (m.ss_obj_value) so that cost → 0 at setpoint.
    """

def default_options():
    """
    Return a dict of Options field overrides to apply as model-level defaults.
    Example: {'nfe_finite': 5, 'sampling_time': 0.2}
    """
```

### Algebraic CV Initial Conditions

`_fix_initial_conditions_at_t0` sets initial conditions at t=0 by fixing all differential state variables directly:

```python
def _fix_initial_conditions_at_t0(block):
    t0 = block.time.first()
    for var in block.state_vars:
        for index in var:
            time_val = index[-1] if isinstance(index, tuple) else index
            if time_val == t0:
                var[index].fix()
```

**Algebraic CVs (e.g., `xD1A`, `xD2B`, `xC`) are NOT fixed at t=0.** Their values at t=0 are implicitly determined by `equations_write` equality constraints — but the framework **deletes all constraint entries at t=0** immediately after `equations_write` runs in `_finite_block_gen` (RADAU: t=0 is not a collocation point). This means algebraic CVs at t=0 have no equality constraints and remain free variables. This is intentional: their actual IC is carried in as the warm-start value from the IC solve and propagated via the unfix→load→refix pattern in `run_MPC.py`.

The steady-state warm-start in `_make_finite_horizon_model` skips `t=0` so that `initialize=...` values are preserved there for non-CV state vars.

The infinite-horizon `endpoint_state_constraints` already handles algebraic CVs generically (`CV[τ=1] == setpoint` for all CVs in `CV_index`, regardless of type).

### IC Consistency Check

After `_fix_initial_conditions_at_t0` fixes all state variables, `_check_ic_consistency(block, tol=1e-4)` (in `tools/initialization_tools.py`) evaluates every active equality constraint at t=0 (skipping `_disc_eq` constraints, which are collocation bookkeeping). If any residual exceeds `tol`, it raises `RuntimeError` listing the top 5 violations by magnitude. This fires automatically at the end of both `_make_finite_horizon_model` and `_make_infinite_horizon_model`.

### IC Builder Functions (`tools/initialization_tools.py`)

The IC data that gets loaded at t=0 is produced by one of two paths, routed by `_build_ic_data(options)`:

**Default path — `_build_ic_data_steady_state(options)`** (`dynamic_initial_conditions=False`):
- Builds the SS NLP model with `fix_slacks=False` (slacks free, so specs can be violated at the start)
- Solves toward `initial_values` targets
- Returns `ScalarData` with consistent values for all variables at t=0
- On IPOPT failure: raises with a readable message explaining that `initial_values` may be infeasible

**Dynamic path — `_build_ic_data_dynamic(options)`** (`dynamic_initial_conditions=True`):
- Builds a model with nfe=1, ncp=1 but deactivates `disc_eq` instead of fixing derivatives to zero — derivatives become free variables determined by the ODE constraints (allowing non-steady-state starting points)
- Applies scalar overrides from `initial_values` to specific state variables
- Fixes all differential state vars at t=0 to their current Pyomo `initialize=` values (set per-index in `variables_initialize`; works for indexed vars like distillation trays)
- Fixes MVs at their `initialize=` values to prevent an underdetermined system
- Solves the resulting feasibility problem for algebraic vars and free derivatives
- On IPOPT failure: calls `_report_constraint_violations` and raises with a clear message

### Slack Variable Squareness Fix

`_make_steady_state_model(m, options, fix_slacks=True)` gains a `fix_slacks` parameter. When `True` (the default, used for the setpoint SS solve), all variables named in `m.slack_index` are fixed to 0 after `variables_initialize`, making the NLP square. The IV path calls it with `fix_slacks=False` so that infeasible-by-spec starting points are still reachable.

### Important: `m.time` is injected by the framework

`variables_initialize` is called with `m.time` already set as a `ContinuousSet`. The user must NOT declare `m.time` in their model. The framework sets:
- Finite block: `m.time = ContinuousSet(bounds=(0, nfe_finite * sampling_time))`
- Infinite block: `m.time = ContinuousSet(bounds=(0, 1))` (transformed domain)
- Plant: `m.time = ContinuousSet(bounds=(0, sampling_time))` (nfe=1)
- Steady-state model: `m.time = ContinuousSet(bounds=(0, 1))` (nfe=1, ncp=1)

---

## Closed-Loop Simulation Flow (`run_MPC.py` → `mpc_loop`)

```
mpc_loop(options)
  │
  ├─ [if initialization_assist] call _assist_initialization_{infinite|finite}(options)
  │   returns a pre-warmed controller
  ├─ [else] InfiniteHorizonController(options) or FiniteHorizonController(options)
  │
  ├─ _create_run_folder(options, resolved_gamma)  ← creates timestamped folder immediately
  │
  ├─ Plant(options)                   # single-step plant, MVs fixed
  │
  ├─ get initial CV values from plant
  │
  └─ for i in range(num_horizons):

        controller.solve()            # IPOPT solve, records cpu time
        [if debug_flag] _report_constraint_violations(controller._model, ...)

        [if lyap_flag]:
          V_current = phi_track[τ=1] on infinite (or finite) block
          first_stage_cost = quadratic cost at t_first_fe on finite block
          update lyap_block.V_prev and lyap_block.first_stage_cost_prev

        # Warm-start plant from controller's full first-FE solution.
        # load_data SKIPS FIXED VARIABLES — MVs and slacks are fixed, so they must
        # be temporarily unfixed to receive the controller's optimal values.
        unfix all plant MVs
        unfix all plant slacks (if slack_index exists)
        plant.interface.load_data(controller.interface.get_data_at_time(new_data_time),
                                  time_points=new_data_time)
        re-fix all plant MVs
        re-fix all plant slacks

        plant.solve()                 # integrate one sampling interval

        extract CV and MV values from plant → io_data_array
        [if safe_run] _save_io_csv(...)

        concatenate plant trajectory onto sim_data

        # Build tf_data: state vars only at sampling_time (ScalarData)
        tf_data = plant state vars at sampling_time

        # Update plant IC for next step (unfix→load→refix at t=0).
        # State vars at t=0 are fixed; load_data skips fixed vars without this pattern.
        plant.interface.load_data(tf_data)   # warm-start t>0 collocation points
        unfix all plant state vars at t=0
        plant.interface.load_data(tf_data, time_points=t0_plant)
        re-fix all plant state vars at t=0

        # Shift controller horizon and load exact plant endpoint as new IC.
        controller.interface.shift_values_by_time(sampling_time)
        unfix all controller state vars at t0_controller
        controller.interface.load_data(tf_data, time_points=t0_controller)
        re-fix all controller state vars at t0_controller
        # Reset extrapolated tail to steady-state guess (avoids Lyapunov infeasibility)
        controller.interface.load_data(_ss_warm_data, time_points=t_last_controller)

  _handle_mpc_results(sim_data, time_series, io_data_array, plant, cpu_time, options,
                      folder_path=run_folder)
```

**Critical**: `DynamicModelInterface.load_data` silently **skips fixed variables**. Every place where data must be loaded into fixed-at-t=0 state vars, or into fixed MVs/slacks, requires the explicit unfix→load→refix pattern shown above. Omitting it means those variables never receive updated values, which causes incorrect initial conditions and plant infeasibility.

---

## Model Construction (`make_model.py`)

This is the most complex file. It builds Pyomo models in stages.

### Infinite-Horizon Model (`_make_infinite_horizon_model`)

1. **Steady-state solve**: Build `_make_steady_state_model` (nfe=1, ncp=1, derivatives fixed to 0). Solve to find `steady_state_values` (setpoints for CVs and MVs) and `ss_obj_value`.
2. **Initial condition solve**: Build another steady-state model, solve toward `initial_values` from model.
3. **Build `m.finite_block`** via `_finite_block_gen()`.
4. **Build `m.infinite_block`** via `_infinite_block_gen()`.
5. **Link blocks** via `_link_blocks()` — state equality constraints at the interface.
6. Create `m.interface = DynamicModelInterface(m.finite_block, m.finite_block.time, clean_model=True)`.
7. Load initial data into finite block; fix state vars at `t=0`.

### Finite-Horizon Block (`_finite_block_gen`)

- Time: `[0, nfe_finite * sampling_time]`
- Collocation: `LAGRANGE-RADAU` with `clean_model='delete'`
- If `ncp_finite > 1`: MVs are reduced to 1 collocation point (piecewise-constant)
- **t=0 constraint deletion**: After calling `m.equations_write(m)`, the block immediately deletes all constraint entries indexed at t=0:
  ```python
  for _con in m.component_objects(pyo.Constraint):
      for _idx in [k for k in _con.keys() if (k[-1] if isinstance(k, tuple) else k) == t0]:
          del _con[_idx]
  ```
  This is necessary because `equations_write` runs post-discretization (when `m.time` is fully populated) and inadvertently creates entries at t=0, which is NOT a RADAU collocation point. Constraints at t=0 would couple algebraic CVs to state variables at that time point — but the state-var ICs are loaded by `run_MPC.py` via unfix→load→refix, not by constraint inference. Leaving these constraints active causes stale algebraic CV values (never updated by `shift_values_by_time`) to corrupt state var ICs after iteration 0.
- If `lyap_flag`: adds `phi_track` as a scalar `pyo.Expression` = sum of quadratic tracking costs at each FE right-endpoint (RADAU collocation points), skipping `t=0`. This is the Riemann-sum contribution of the finite block to the total Lyapunov measure.
  - For finite-only: also adds `V_prev`, `first_stage_cost_prev`, and `lyap_stability_constraint` on the block.
  - For infinite-horizon: the Lyapunov constraint is placed on the parent `ConcreteModel` by `_make_infinite_horizon_model` (see below).
- `m.state_vars`, `m.deriv_vars` are stored on the block

### Infinite-Horizon Block (`_infinite_block_gen`)

- Time: `[0, 1]` (transformed domain, maps to `[sampling_time, ∞)` physically)
- Collocation: `LAGRANGE-LEGENDRE` with `clean_model='delete'`
- **Gamma computation**: If `options.gamma is None`, calls `_compute_gamma_from_collocation(nfe_infinite, ncp_infinite, sampling_time)` → `atanh(tau_1) / sampling_time` where `tau_1` is the first collocation point. This ensures `tanh(gamma * sampling_time) == tau_1`.
- **Added variables/constraints**:
  - `phi`: primary terminal cost integral (ODE `dphidt`)
  - `phi[0].fix(0)`; `terminal_cost` constraint defines the ODE dynamics
  - `phi_track`: quadratic tracking Lyapunov integral (ODE `dphidt_track`)
  - `phi_track[0].fix(0)`; `phi_track_ode` defines its dynamics (starts at 0, independent of finite block)
  - Terminal constraints (when `terminal_constraint_type != 'none'`): `diff_terminal_constraints` for differential CVs (direct τ=1 access); `alg_terminal_constraints` for algebraic CVs/MVs (Lagrange extrapolation to τ=1). For `'soft'`, an `infinite_terminal_soft_penalty` Expression is built instead and added to the objective by the controller.
- **Derivative transformation** (`_transform_model_derivatives`): Replaces all `dxdt[t]` in user constraints with `(gamma/sampling_time * (1 - t²)) * dxdt[t]`. Skips `*_disc_eq`, `terminal_cost`, and `phi_track_ode` constraints.

### Infinite-Horizon Lyapunov Constraint (top-level model, set in `_make_infinite_horizon_model`)

When `lyap_flag=True` and `infinite_horizon=True`, the combined Lyapunov stability constraint is placed on the parent `ConcreteModel` (not on either block):
```python
m.lyap_stability_constraint:
    finite_block.phi_track + infinite_block.phi_track[τ=1] - V_prev
    <= -lyap_delta * first_stage_cost_prev
```
`V_prev` and `first_stage_cost_prev` live on the top-level model `m`.
`finite_block.phi_track` is the Riemann-sum Expression; `infinite_block.phi_track[τ=1]` is the ODE integral.

The constraint form is controlled by `lyap_constraint_type`:
- `'hard'` — strict inequality constraint as shown above.
- `'soft'` — adds `m.lyap_slack` (`NonNegativeReals` Var) to the RHS, allowing the constraint to be violated at cost `lyap_soft_weight * lyap_slack` in the objective. Useful when a short finite horizon makes the hard constraint locally infeasible after a shift.
- `'none'` — `phi_track` and `V_prev` infrastructure is still built, but no constraint is added.

### Time Transformation (key math)

Physical ODE: `dx/dt = f(x, u)`
Transformed domain `τ = tanh(γt) ∈ [0,1)`:
`dx/dτ = (sampling_time/γ) * (1/(1-τ²)) * dx/dt`

So in user constraints, `dxdt[τ]` gets replaced by the scaled version. The terminal cost and phi_track ODEs explicitly include the Jacobian factor `(γ/sampling_time * (1-τ²))`.

### Linking Blocks (`_link_blocks`)

For each state variable: `finite_block.var[t_end] == infinite_block.var[t_start]`

---

## Controller Classes (`controllers.py`)

### `Controller` (base)
- Holds `self._model` (Pyomo ConcreteModel), `self._solver` (cold IPOPT), `self._warm_solver` (warm-start IPOPT), `self._initialized = False`
- **Cold solver** (`_ipopt_solver`): used only on the first `solve()` call. Standard IPOPT settings with `bound_relax_factor=0`.
- **Warm solver** (`_ipopt_warm_solver`): used on all subsequent `solve()` calls after shift+load. Same settings as cold but with `bound_push=1e-8` and `bound_frac=1e-8` to keep the warm-started solution close to bounds. Note: `warm_start_init_point` is intentionally omitted because Pyomo has no Suffix setup for dual multipliers — enabling it would pass zero duals to IPOPT, which is worse than the default.
- `solve()`: dispatches to cold or warm solver via `_initialized` flag, checks termination condition, records `self.last_solve_time` (CPU via `resource.getrusage`), then sets `_initialized = True`
- `__getattr__`: falls through to `self._model` — so `controller.finite_block`, `controller.interface`, etc. all work

### `InfiniteHorizonController(options, data=None)`
1. Calls `_make_infinite_horizon_model`
2. Sets `m.finite_block.stage_cost_index = CV_index + MV_index`
3. Builds objective (one of three modes):
   - **custom_objective**: `Σ_t [stage_cost(t) - ss_obj_value]` over finite elements + `beta/sampling_time * phi[τ=1]`
   - **terminal_cost_riemann**: quadratic stage cost + Riemann sum over infinite block
   - **standard quadratic**: quadratic stage cost + `beta/sampling_time * phi[τ=1]`
4. Dumps model to `model_output.txt` (always, at construction time)
5. Performs initial solve

### `FiniteHorizonController(options, data=None)`
- Same pattern but single block, no infinite block
- Objective sums quadratic (or custom) stage cost over finite elements only
- Also dumps model to `model_output.txt` at construction time

---

## Plant Class (`plant.py`)

```python
Plant(options)
```
- Internally copies options with `nfe_finite=1, infinite_horizon=False`
- Builds a `_make_finite_horizon_model` with 1 finite element
- Fixes all MVs (`var.fix()`)
- If `ncp_finite > 1`: deactivates `{mv_name}_interpolation_constraints`
- Adds dummy objective `expr=1`
- `solve()`: IPOPT solve; on failure, calls `_report_constraint_violations` if `debug_flag=True`
- `__getattr__`: falls through to `self._model`

---

## Model Validator (`check_model.py`)

A standalone validator for user `model.py` files. Runs three levels of checks:

**Level 1 — Static AST:**
- Valid Python syntax
- `variables_initialize()` and `equations_write()` are defined
- `custom_objective()` and `default_options()` detected (WARN if absent)
- `variables_initialize()` does not declare `m.time` (framework injects it)
- Both required functions return `m`

**Level 2 — Import:**
- Module loads without error

**Level 3 — Dynamic Pyomo:**
- `variables_initialize()` runs on a real `ConcreteModel` with `m.time` set
- `m.MV_index`, `m.CV_index`, `m.setpoints`, `m.initial_values` are declared
- At least one `DerivativeVar` is present
- `m.MV_index` and `m.CV_index` are non-empty
- `equations_write()` runs after collocation discretization

**Usage:**
```bash
check-model path/to/model.py            # full check (static + import + dynamic)
check-model path/to/model.py --no-dynamic  # skip Pyomo dynamic check
```

**Programmatic API:**
```python
from infeNMPC.check_model import validate
report = validate(path, dynamic=True)
print(report.ok)       # True if no failures
print(report.failed)   # number of FAIL results
```

---

## Debug Tools (`tools/debug_tools.py`)

```python
_report_constraint_violations(model, n=10, tol=1e-6, label="")
```

Traverses all active constraints (including nested blocks) and prints the `n` most violated. Used automatically when `options.debug_flag=True`:
- After every successful controller solve (in `run_MPC.py`)
- After a failed controller solve (in `run_MPC.py`, before re-raising)
- After a failed plant solve (in `plant.py`, before re-raising)

---

## Data Flow: Controller ↔ Plant

The `DynamicModelInterface` (from `pyomo.contrib.mpc`) is the bridge. Both `controller.interface` and `plant.interface` point to `DynamicModelInterface` instances wrapping their respective Pyomo models.

Key operations in each MPC iteration:
```python
# 1. Get controller's optimal MV at sampling_time
ts_data = controller.interface.get_data_at_time(options.sampling_time)
input_data = ts_data.extract_variables([controller.finite_block.Fa0, ...])

# 2. Load MVs into plant
plant.interface.load_data(input_data, time_points=new_data_time)

# 3. Solve plant
plant.solve()

# 4. Get plant state at end of interval (state vars only, not MVs)
tf_data = ScalarData({key: value, ...})

# 5. Update plant IC for next step
plant.interface.load_data(tf_data)

# 6. Shift controller horizon and update IC
controller.interface.shift_values_by_time(options.sampling_time)
controller.interface.load_data(tf_data, time_points=t0_controller)
```

Variable keys have the form `"Ca[*]"` (scalar-in-time) or `"x[1,*]"` (spatially indexed). These are generated by `_get_variable_key_for_data(model, name)` in `tools/indexing_tools.py`.

---

## Lyapunov Stability (`lyap_flag=True`)

**Finite-horizon only** (`infinite_horizon=False`): `V_prev`, `first_stage_cost_prev`, and `lyap_stability_constraint` live on the finite block (= `controller._model`):
```python
m.lyap_stability_constraint:
    phi_track - V_prev <= -lyap_delta * first_stage_cost_prev
```
where `phi_track` is a scalar `pyo.Expression` = Riemann sum of quadratic tracking costs at FE endpoints.

**Infinite-horizon** (`infinite_horizon=True`): `V_prev`, `first_stage_cost_prev`, and `lyap_stability_constraint` live on the top-level `ConcreteModel` (= `controller._model`):
```python
m.lyap_stability_constraint:
    finite_block.phi_track + infinite_block.phi_track[τ=1] - V_prev
    <= -lyap_delta * first_stage_cost_prev
```
The total Lyapunov measure combines the Riemann sum over the finite block and the ODE integral over the infinite block.

`phi_track` (the quadratic tracking measure) always uses `Σ c[i] * (var[i] - setpoint[i])²` with `stage_cost_weights` — independent of `custom_objective`.

At each MPC iteration (in `run_MPC.py`):
```python
# lyap_block = controller (proxies to controller._model) in both cases
if infinite_horizon:
    V_current = pyo.value(controller.finite_block.phi_track) \
              + pyo.value(controller.infinite_block.phi_track[tau_last])
else:
    V_current = pyo.value(controller.phi_track)
first_stage_cost = Σ c[j] * (var[j][t_first_fe] - setpoint[j])²  # at first FE endpoint
lyap_block.V_prev.set_value(V_current)
lyap_block.first_stage_cost_prev.set_value(first_stage_cost)
```

`lyap_delta = 0.01` is a loose constraint (1% decrease required per iteration). `lyap_delta = 0.99` is tight (99% decrease required per iteration). Constraint form: `V_k - V_{k-1} ≤ -δ * L_{k-1}`.

---

## Results Output (`data_save_and_plot.py`)

Results are saved to a timestamped folder created once at the start of `mpc_loop`:
```
Results/
  {model_name}/
    {YYYY-MM-DD}/
      {HH-MM-SS}/
        run_config.txt    ← all Options fields + resolved gamma + git branch
        io_data.csv       ← CV and MV trajectories (Time, Ca, Cb, Fa0, ...)
        sim_data.csv      ← full DAE trajectory from plant
        plot_results.py   ← auto-generated, standalone plotting script
        CV_1_Ca.png, CV_2_Cb.png, MV_1_Fa0.png  ← trajectory plots
```

`model_name` comes from `options.model_name` (set to the module's `__name__` last component by `Options.for_model_module`). Special characters are replaced with `_`.

`run_config.txt` is written immediately when the folder is created (before the loop starts). The `gamma` field shows the resolved value even when `options.gamma=None` (auto-computed).

For `safe_run=True`, `io_data.csv` is overwritten after every iteration into the same timestamped folder.

### Auto-generated `plot_results.py`

`_write_plot_script()` is called at the end of every run and deposits a self-contained Python script into the results folder. It is pre-populated with:
- Active `PLOTS` entries for every CV and MV (using their `CV_display_names` / `MV_display_names` as LaTeX labels)
- Commented-out entries for all remaining `sim_data.csv` columns (easy to uncomment for state-variable diagnostics)
- `SETPOINTS` values from `plant.steady_state_values` for each CV
- Both PDF (vector) and PNG (raster) output per figure

Run it standalone: `python Results/.../plot_results.py`

### Key functions

| Function | Purpose |
|---|---|
| `_create_run_folder(options, resolved_gamma)` | Create `Results/model/date/time/`, write `run_config.txt`, return path |
| `_handle_mpc_results(..., folder_path)` | Save CSVs + figures + `plot_results.py` into `folder_path` |
| `_write_plot_script(folder_path, options, plant, sim_data)` | Generate standalone plotting script |
| `_get_results_folder(options)` | Legacy: returns `Results/{model_name}/` (no timestamp) |

---

## Tools

### `collocation_tools._compute_gamma_from_collocation(nfe, ncp, sampling_time)`
Builds a throwaway Pyomo model, discretizes it with LAGRANGE-LEGENDRE, reads the first non-boundary collocation point `tau_1`, returns `atanh(tau_1) / sampling_time`.

### `indexing_tools`
- `_get_variable_key_for_data(model, name)` → `"Ca[*]"` or `"x[1,2,*]"`
- `_parse_indexed_name("x[1,2]")` → `("x", (1, 2))`
- `_add_time_indexed_expression(model, var_name, t)` → Pyomo expression `m.x[1, t]`
- `_get_derivative_and_state_vars(model)` → `(set of DerivativeVars, set of state Vars)`
- `_get_disc_eq_time_points(m)` → sorted list of collocation time points

### `debug_tools`
- `_report_constraint_violations(model, n=10, tol=1e-6, label="")` → print top-N most violated constraints

### `initialization_tools`
Progressive warm-start: solve sequence of controllers starting at `initialization_assist_sampling_time_start` (e.g., 10.0), reducing by factor 0.9 until reaching `options.sampling_time`. Each solve's solution seeds the next. Used for stiff systems that won't converge from cold-start at the target sampling time.

---

## IPOPT Configuration

All solvers are configured identically via `_ipopt_solver()`:
```python
solver.options['linear_solver'] = 'ma57'
solver.options['OF_ma57_automatic_scaling'] = 'yes'
solver.options['max_iter'] = 500
solver.options['halt_on_ampl_error'] = 'yes'
solver.options['bound_relax_factor'] = 0
```

---

## Pyomo/IDAES Conventions

- **`clean_model='delete'`**: Passed to `dae.collocation.apply_to()`. Removes non-collocation variable/constraint entries after discretization. Requires forked Pyomo (`oag00002/pyomo`, branch `mpc-rewrite`).
- **`DynamicModelInterface(..., clean_model=True)`**: Uses sparse data-access paths for compatibility with `clean_model='delete'`.
- **Discretization scheme**:
  - Finite block → `LAGRANGE-RADAU`
  - Infinite block → `LAGRANGE-LEGENDRE`
  - Steady-state → `LAGRANGE-RADAU` (nfe=1, ncp=1)
- **MV collocation reduction**: When `ncp_finite > 1`, MVs are reduced to `ncp=1` (piecewise-constant) via `discretizer.reduce_collocation_points(...)`.

---

## Example: CSTR with Lyapunov Constraint (canonical working example)

**System**: Isothermal CSTR, A→B consecutive reaction
**Location**: `examples/lyap_flag_examples/enmpc_cstr/`

```
small_cstr_model.py:
  MV_index = ["Fa0"]          (feed rate of A, bounds [10, 20] kmol/h)
  CV_index = ["Ca", "Cb"]     (concentrations in mol/L)
  ODEs: dCa/dt = Fa0/V*(Caf - Ca) - k*Ca
        dCb/dt = Fa0/V*(0 - Cb) + k*Ca
  custom_objective: lambda m, t: -m.Fa0[t]*(2*m.Cb[t] - 0.5)
  setpoints: Ca=0.5, Cb=0.5
  initial_values: Ca=0, Cb=0

run.py:
  options = Options.for_model_module(small_cstr_model,
    num_horizons=100, sampling_time=1,
    nfe_finite=2, ncp_finite=1,
    nfe_infinite=5, ncp_infinite=1,
    infinite_horizon=True, terminal_constraint_type='none',
    custom_objective=True,
    stage_cost_weights=[1, 1, 1/600], beta=1,
    lyap_flag=True, lyap_delta=0.01,
    save_data=True, save_figure=True,
  )
  mpc_loop(options)
```

Run from repo root: `python examples/lyap_flag_examples/enmpc_cstr/run.py`

---

## Example: Two-Column Ternary Distillation (new large-scale example)

**System**: Two distillation columns in series separating A/B/C
**Locations**:
- `examples/standard/ternary_distillation_double_column/` — standard (no Lyapunov)
- `examples/lyap_flag_examples/enmpc_ternary_distillation/` — Lyapunov support (`lyap_flag=False` by default)

```
ternary_distillation_model.py:
  MV_index = ["VB1", "LT1", "D1", "B1", "VB2", "LT2", "D2", "B2"]
  CV_index = ["xD1A", "xD2B", "xC"]   (algebraic CVs, not differential)
  246 state variables: x1[tray,comp], M1[tray] (Col 1); x2[tray,comp], M2[tray] (Col 2)
  Francis Weir formula liquid dynamics
  custom_objective: minimize feed + energy cost minus product revenue
  slack_index: ["xD1A_eps", "xD2B_eps", "xC_eps", "M1_eps", "M2_eps"]

run.py (lyap_flag_examples version, confirmed working as of 2026-04-15):
  options = Options.for_model_module(ternary_distillation_model,
    num_horizons=5,         # set to 5 for testing; restore to 10+ for production
    sampling_time=1,
    nfe_finite=1, ncp_finite=3,
    infinite_horizon=False,
    custom_objective=True,
    lyap_flag=True, lyap_delta=0.01,
    lyap_constraint_type='soft', lyap_soft_weight=1.0,
    terminal_constraint_type='soft',
    debug_flag=True, safe_run=True, tee_flag=True,
  )
```

This model uses **algebraic CVs** (`xD1A`, `xD2B`, `xC` are computed from tray compositions via equality constraints). They are NOT fixed at t=0 by `_fix_initial_conditions_at_t0` — only differential state vars are fixed. The t=0 constraint deletion in `_finite_block_gen` removes the `xD1A_def[0]`, `xD2B_def[0]`, `xC_def[0]` entries that would otherwise couple algebraic CVs to state vars at t=0. It also demonstrates `m.slack_index` — an optional set for slack variables used for constraint softening.

**Plant options override**: `plant_options` in `Plant.__init__` sets `lyap_flag=False` to eliminate the `lyap_slack` DOF from the plant model (plant only needs to integrate the ODE; Lyapunov slack is a controller artifact).

---

## Known Issues / State at Last Review (2026-04-15)

- Branch: `shifting-behavior`
- Most examples are not confirmed working; the CSTR lyap_flag example is the primary test case
- **Ternary distillation lyap example confirmed working** (5 iterations, all IPOPT-optimal, constraint violations at 1e-12 to 1e-15 numerical precision). See LLMJOURNAL.md "Plant Squareness & Warm-Start" entry for full history.
- `model_output.txt` is written to the working directory at controller construction time (`controllers.py`). This is debug behavior.
- `gamma` in the results folder path reflects `options.gamma` (which stays `None` until the block stores it on `infinite_block.gamma` — the folder path may show `gamma_None` even after auto-computation)
- **Pending cleanup items** (not yet addressed):
  - `plant.py`: unconditional `pprint` to `plant_model.txt` in `Plant.__init__` should be gated by `options.debug_flag` or removed
  - `run_MPC.py`: pprint dumps for iterations 0 and 1 (controller post-solve and plant pre-solve) gated by `debug_flag and i <= 1` — acceptable for debugging, but can be removed for production
  - `examples/lyap_flag_examples/enmpc_ternary_distillation/run.py`: `num_horizons=5` was set for testing; should be restored to 10+ for production runs
  - **Infinite block (`_infinite_block_gen`)**: Uses LAGRANGE-LEGENDRE where both `t=0` AND `t=t_last` are non-collocation points. The same `equations_write` spurious-constraint issue almost certainly exists there. Not yet investigated.
- **IC framework added (2026-04-09):** `_make_steady_state_model` has `fix_slacks=True` param; `_check_ic_consistency`, `_build_ic_data_steady_state`, `_build_ic_data_dynamic`, `_build_ic_data` added to `tools/initialization_tools.py`; `dynamic_initial_conditions` option added. See LLMJOURNAL.md for full details.
