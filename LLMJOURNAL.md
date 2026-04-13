# Next Tasks

## Terminal Constraint Overhaul (2026-04-10)

**Branch:** `fix-terminal-constraints-2`

Replaced the old `endpoint_constraints: bool` option with a proper terminal constraint framework
covering all four combinations of horizon type × constraint style.

---

### Motivation

The old `endpoint_constraints=True` had several correctness and completeness gaps:

1. **Wrong physics for algebraic CVs/MVs.** With LAGRANGE-LEGENDRE collocation, τ=1 is not a
   quadrature point. Differential state variables are pinned at τ=1 by Pyomo's internal
   continuity equations. Algebraic variables (CVs without a `DerivativeVar`, and all MVs) have
   no such equation — their τ=1 index is a free variable with no governing constraint. The old
   code added `xD1A[1] == ss_val` directly, which constrained a free variable that had no
   connection to the actual physics at the collocation points.
2. **MVs were never constrained** at the terminal time.
3. **No soft-constraint option** (penalty-only alternative).
4. **No terminal constraints in finite-horizon controllers** at all.

---

### New Options (`infNMPC_options.py`)

`endpoint_constraints: bool` removed. Three new fields added:

| Option | Default | Meaning |
|--------|---------|---------|
| `terminal_constraint_type` | `'hard'` | `'hard'` = equality; `'soft'` = penalty; `'none'` = off |
| `terminal_constraint_variables` | `'cvmv'` | `'cv'`, `'mv'`, or `'cvmv'` |
| `terminal_soft_weight` | `1.0` | Scalar multiplier on soft penalty (per-var weights from `stage_cost_weights`) |

All ~10 example `run.py` files updated: `endpoint_constraints=True` → `terminal_constraint_type='hard'`,
`endpoint_constraints=False` → `terminal_constraint_type='none'`.

---

### Polynomial Extrapolation for Algebraic Variables

For the infinite-horizon block, algebraic CVs and all MVs require Lagrange interpolation to
compute the terminal value:

```
var_terminal ≈ Σ_j  L_j(1) · var[τ_j]
```

where `τ_j` are the collocation points in the last finite element, and `L_j(1)` are the Lagrange
basis polynomials evaluated at 1 (pure numerical scalars, computed once).

Two helpers added to `tools/collocation_tools.py`:

- `_get_last_fe_pts(block)` — returns sorted collocation points strictly inside the last FE
- `_lagrange_coeffs_at_endpoint(colloc_pts, t_end)` — computes `L_j(t_end)` for each point

---

### Implementation Details

**`make_model.py`** — two module-level helpers added:
- `_terminal_vars(m, options)` — builds the list of variable names to constrain based on
  `terminal_constraint_variables`
- `_terminal_weights(options, var_names, stage_cost_index)` — maps each variable to its weight
  from `stage_cost_weights` by position in `stage_cost_index`

**`_finite_block_gen`** — terminal constraint block added after `equations_write`, guarded by
`not options.infinite_horizon` (so it fires only for pure finite-horizon controllers, not when
the finite block is embedded inside an infinite-horizon model). RADAU collocation includes
the right FE endpoint, so `m.time.last()` is a valid collocation point for all variables.
- `'hard'` → `m.finite_terminal_constraints` (one equality per variable)
- `'soft'` → `m.finite_terminal_soft_penalty` Expression (quadratic penalty sum)

**`_infinite_block_gen`** — `endpoint_state_constraints` block replaced. Variables are split into:
- **Differential CVs** (`var in state_vars`) → `m.diff_terminal_constraints`: direct equality
  `CV[τ=1] == ss_val` (Pyomo continuity equations already pin this value)
- **Algebraic CVs / MVs** (`var not in state_vars`) → `m.alg_terminal_constraints`: equality
  on the Lagrange extrapolation expression
- For `'soft'`: `m.infinite_terminal_soft_penalty` Expression built from the same extrapolation
  logic; no equality constraints added

**`controllers.py`** — all three infinite-horizon objective branches and both finite-horizon
branches add `*_terminal_soft_penalty * terminal_soft_weight` to their return value when
`terminal_constraint_type='soft'`.

## ~~Fix phi_track in the finite horizon block~~ (Done)

**Done.** `phi_track` in the finite block is now a scalar `pyo.Expression` (Riemann sum over
FE endpoints). For infinite-horizon, the Lyapunov constraint on the top-level `ConcreteModel`
combines both blocks: `finite_block.phi_track + infinite_block.phi_track[τ=1]`.
Results folder restructured to `Results/{model}/{date}/{time}/` with `run_config.txt`.

---

## Fix phi_track in the finite horizon block (original description)

**File:** `infeNMPC/make_model.py` — `_finite_block_gen`

**Problem:** `phi_track` is currently built as a continuous integral via ODE + DerivativeVar,
discretized over all collocation points (Gaussian quadrature). It should mirror the objective:
a simple sum of the quadratic tracking cost evaluated at each finite element endpoint.

**What to change:**
- Remove `phi_track` as a `pyo.Var` with `DerivativeVar` and `phi_track_ode` constraint.
- Replace with a scalar `pyo.Var` (or `pyo.Expression`) equal to the sum:
  ```
  phi_track_val = sum(tracking_cost(t) for t in get_finite_elements(), t != t_first)
  ```
- Update the Lyapunov stability constraint (`lyap_stability_constraint`) to reference the new scalar.
- Update `run_MPC.py` where `V_current = pyo.value(lyap_block.phi_track[lyap_block.time.last()])` —
  this needs to read the scalar instead of an indexed Var.

**MV note:** With LAGRANGE-RADAU, FE right-boundary endpoints are Radau collocation points
and MVs exist there (via `reduce_collocation_points` interpolation constraints), so summing
over `get_finite_elements()` is safe. If the scheme ever changes to LEGENDRE, use the first
collocation point of each element instead.

---

## Tune lyap_delta for ternary distillation

~~The ternary distillation `run.py` currently has `lyap_delta=0.99`, which is far too tight.~~
**Done** — changed to `lyap_delta=0.01`. Tighten further once results are reviewed.

---

## Fix initial condition solve — always solve the iv model

**Files:** `infeNMPC/make_model.py` — `_make_finite_horizon_model` and `_make_infinite_horizon_model`

### The bug

Both functions contain this shortcut:

```python
initial_value_vars = list(m_iv.initial_values.index_set())
if all(var in m_iv.CV_index for var in initial_value_vars):
    initial_data = ScalarData({
        _get_variable_key_for_data(m_iv, cv): pyo.value(m_iv.initial_values[cv])
        for cv in initial_value_vars
    })
else:
    ...
    initial_data = _solve_steady_state_model(m_iv, m_iv_target, options, label="iv")
```

When `initial_values` is indexed only by `CV_index` (as in every current example), the
iv solve is skipped entirely. The resulting `initial_data` contains **only the CV values** —
nothing for the non-CV state variables.

### Why this is wrong

`m.interface.load_data(initial_data, time_points=t_first)` then loads this sparse data into
the finite block at `t=0`. State variables not in `initial_data` are left at their Pyomo
`initialize=...` defaults. Those defaults are model-author-supplied guesses — **not**
generally consistent with the algebraic constraints or ODEs at `t=0`.

For the ternary distillation model, the CVs (`xD1A`, `xD2B`, `xC`) are algebraic — they are
not themselves differential state variables. The 246 actual state variables (`x1[tray,comp]`,
`M1[tray]`, `x2[tray,comp]`, `M2[tray]`) all have `initialize=...` values that may not be
consistent with those CV targets. `_fix_initial_conditions_at_t0` correctly handles the one
tray entry per CV that is linked through a defining constraint (e.g. `x1[41,1,0]` through
`xD1A_def`), but the remaining ~240 state variables at `t=0` are just fixed at their guesses.

There is also no good reason for the shortcut in the first place. Solving the iv model when
`initial_values` only contains CVs is harmless: the solver targets those CV values and returns
a fully consistent `ScalarData` for all variables, including all the differential states. That
gives every state variable a physically consistent IC instead of an independent guess.

### Fix

Remove the `if all(var in m_iv.CV_index ...)` branch entirely. Always build `m_iv_target`
from `m_iv.initial_values` and call `_solve_steady_state_model` regardless of whether the
initial variables are CVs or not. The returned `initial_data` will then contain consistent
values for all variables, which `load_data` can spread across `t=0`.

```python
m_iv_target = ScalarData({
    _get_variable_key_for_data(m_iv, var): pyo.value(m_iv.initial_values[var])
    for var in initial_value_vars
})
initial_data = _solve_steady_state_model(m_iv, m_iv_target, options, label="iv")
```

Apply this change in both `_make_finite_horizon_model` and `_make_infinite_horizon_model`.

---

## IC Consistency & Dynamic Initial Conditions Framework (2026-04-09)

**Branch:** `main`

Four related bugs/gaps in the initial condition (IC) machinery were addressed in a single session. All changes preserve backward compatibility — existing examples run unchanged with no options modifications.

---

### Problem 1: Slack variables made the SS NLP non-square

Models with `m.slack_index` (e.g., ternary distillation declares `xD1A_eps`, `xD2B_eps`, `xC_eps`, `M1_eps`, `M2_eps`) have unbounded non-negative slack variables that, when left free in the steady-state model, cause the NLP to be under-determined (more free variables than constraints). IPOPT would either hang or return spurious solutions.

**Fix:** `_make_steady_state_model` in `make_model.py` gained a `fix_slacks=True` parameter. When True, after `variables_initialize` and before discretization, all variables named in `m.slack_index` are fixed to 0. The setpoint SS solve uses `fix_slacks=True` (the default) so the true economic optimum is found with all specs enforced. The IV solve uses `fix_slacks=False` so the initial point can be off-spec (slacks absorb constraint violations) — important when the user's starting conditions don't yet satisfy product purity specs.

---

### Problem 2: No check that ICs are actually consistent

After loading `initial_data` at t=0 and fixing state variables, nothing verified that the algebraic constraints at t=0 were satisfied. If the IC solve returned a slightly infeasible point (or if the user's `initialize=` values were inconsistent with the equations), the error would only appear later as a cryptic IPOPT convergence failure in the controller.

**Fix:** New function `_check_ic_consistency(block, tol=1e-4)` in `tools/initialization_tools.py`. It iterates all active equality constraints on the block, evaluates each one at t=0 (skipping `_disc_eq` constraints which are collocation bookkeeping rather than physics), and collects any whose residual exceeds `tol`. If violations are found, it raises `RuntimeError` listing the top 5 by magnitude — giving the user a direct pointer to which physical equations are violated at the initial time. Called at the end of both `_make_finite_horizon_model` and `_make_infinite_horizon_model`, after `_fix_initial_conditions_at_t0`.

---

### Problem 3: No dynamic (non-steady-state) IC path

The only IC path was to solve a steady-state NLP toward `initial_values` — derivatives are fixed to zero, the solver finds the nearest consistent static point. This is the right approach most of the time, but it prevents users from starting at a point that is not near any steady state (e.g., an empty column, a CSTR at zero concentration, or an explicitly prescribed transient initial state).

**Fix:** New `dynamic_initial_conditions: bool = False` option (added to `Options` in `infNMPC_options.py`). When True, `_build_ic_data_dynamic` (in `tools/initialization_tools.py`) is called instead of the SS-NLP path. It:

1. Builds the same nfe=1, ncp=1 model structure
2. Applies scalar overrides from `initial_values` to state var initial values (the `initialize=` Pyomo parameter already encodes per-index ICs for indexed vars like distillation trays)
3. Discretizes (nfe=1, ncp=1, LAGRANGE-RADAU)
4. Calls `equations_write` so all algebraic constraints are present
5. **Deactivates `disc_eq`** rather than fixing derivatives to zero — this leaves derivatives as free variables, their values to be determined by the ODE constraints. This is the key distinction: derivatives can be non-zero at the starting point.
6. Fixes all differential state variables at t_first to their current (Pyomo `initialize=`) values
7. Fixes MVs at their current values to prevent an underdetermined system
8. Sets a zero objective (feasibility problem)
9. Solves with IPOPT; on failure, calls `_report_constraint_violations` and raises with a clear message

This path supports both scalar models (CSTR: Ca, Cb) and spatially-indexed models (distillation: x1[tray, comp]). For indexed vars, the per-tray `initialize=` values already encode the desired IC — no enumeration in `initial_values` is required. `initial_values` only provides optional scalar overrides (e.g., a CV target like `xD1A=0.90`).

---

### Problem 4: Generic error message on IV solve failure

When the steady-state IV NLP failed, the user received a raw IPOPT termination message with no context about what went wrong or why.

**Fix:** `_build_ic_data_steady_state` (the refactored SS-NLP IV path) wraps the `_solve_steady_state_model` call in a try/except that re-raises with a human-readable message: `"Initial condition solve failed — the specified initial_values may define an infeasible or inconsistent starting point. Check that initial_values are achievable under the model equations."` The original exception is chained so the IPOPT details are still available.

---

### Code organization

All new IC logic lives in `tools/initialization_tools.py` to avoid cluttering `make_model.py`. The four new functions are:

| Function | Purpose |
|---|---|
| `_check_ic_consistency(block, tol)` | Post-fix consistency check at t=0 |
| `_build_ic_data_steady_state(options)` | SS-NLP IV path (existing behavior, extracted + improved error) |
| `_build_ic_data_dynamic(options)` | Dynamic (non-SS) IV path |
| `_build_ic_data(options)` | Router: dispatches to one of the above based on `options.dynamic_initial_conditions` |

`make_model.py` changes are minimal: `_make_steady_state_model` gains the `fix_slacks` param; both orchestrators replace their inline IV solve block with `initial_data = _build_ic_data(options)` and append `_check_ic_consistency(...)`.

Circular import avoided: `initialization_tools.py` imports from `make_model` and `model_equations` inside the function bodies (not at module level), since `make_model` already imports from `tools/`.
