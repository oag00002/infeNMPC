# Next Tasks

## Solver Hardening + `revive_run` Flag (2026-04-29)

**Branch:** `binary-distillation`

### Problem: random "restoration phase converged to small primal infeasibility" failures

The MPC loop was failing randomly around 90% of the way through long runs with:
```
Restoration phase converged to a point with small primal infeasibility.
```
IPOPT output showed restoration iterating at constant `inf_pr ≈ 2.38e-7`,
`inf_du ≈ 1e-14`, `mu ≈ 10^{-8.8}`, step halved every iteration (filter rejection),
constant objective. The solution was essentially already found — the failure was purely
numerical.

**Root cause:** `bound_relax_factor = 0` in the warm solver is too strict. After
`shift_values_by_time + IC load`, the warm-started point has a floating-point residual
(~2.38e-7) that lies just outside a bound. Restoration minimises constraint violation
but cannot reach zero because the restoration problem's optimum requires a tiny bound
violation — which `bound_relax_factor = 0` disallows.

`bound_relax_factor = 0` was historically intentional: it prevented IPOPT from returning
slightly out-of-bound values that then caused downstream `set_value()` failures on
`NonNegativeReals` slack variables (e.g., slack set to `-1e-9`). The fix to allow
`bound_relax_factor = 1e-8` safely is to add a post-solve bound-clip step.

---

### Changes

**`infeNMPC/make_model.py`** — `_ipopt_warm_solver`:
- Removed the stale `bound_relax_factor = 0` line (was overriding itself to 0 before new code added 1e-8 below it — cleaned up)
- Added `bound_relax_factor = 1e-8`
- Added `acceptable_tol = 1e-5`, `acceptable_constr_viol_tol = 1e-4`, `acceptable_obj_change_tol = 1e20`, `acceptable_iter = 3` (secondary safety net: if 3 consecutive iterations are all within these looser tolerances, IPOPT exits as "Solved to Acceptable Level" = `TerminationCondition.optimal`)
- Added `_ipopt_revive_solver()`: same as warm solver plus `tol = 1e-6`

**`infeNMPC/tools/solver_tools.py`** — **new file**:
- `_clip_to_bounds(model)`: iterates all Var data objects; clips any value outside `[lb, ub]` back to the bound. Called after every successful controller solve. Handles the downstream concern from `bound_relax_factor != 0`.
- `_attempt_revive(controller, i, options, caught_exc)`: retry logic. If `revive_run > 0` and TC is not `infeasible` or `maxIterations`, retries up to `revive_run` times using `controller.revive_solve()` with clear print messages. Re-raises caught_exc or the last revive failure on exhaustion.
- `_SKIP_REVIVE = frozenset({infeasible, maxIterations})`

**`infeNMPC/controllers.py`**:
- Imports `_ipopt_revive_solver`
- `Controller.__init__`: adds `self._revive_solver = _ipopt_revive_solver()`
- `Controller.solve()` except block: stores `self.last_tc = results.solver.termination_condition` before re-raising (so `_attempt_revive` can inspect it)
- New `Controller.revive_solve()` method: calls `_revive_solver`, sets `self.last_solve_time = revive_time` (assignment, not accumulation — failed initial solve not counted)

**`infeNMPC/infNMPC_options.py`**:
- `revive_run: int = 0` field + docstring

**`infeNMPC/run_MPC.py`**:
- Imports `_clip_to_bounds`, `_attempt_revive` from `tools.solver_tools`
- Replaced the old try/except (which always re-raised) with:
  ```python
  try:
      controller.solve()
  except RuntimeError as exc:
      _attempt_revive(controller, i, options, exc)
  _clip_to_bounds(controller._model)
  cpu_time.append(controller.last_solve_time)
  ```
- `cpu_time.append` moved to after `_clip_to_bounds`

---

## Tracking Objective: Slack Penalty & AMPL-Consistent Weights (2026-04-21)

**Branch:** `fix-tracking`

### Motivation

When `objective='tracking'`, the model's soft-constraint slack variables (`m.slack_index`) had
no penalty in the objective — they were free non-negative variables with zero cost. IPOPT could
drive them to any value to satisfy the soft constraints without any resistance. In the AMPL
reference (`double_col_dyn.run`, `column2_dynamic_soft.mod`), all slack variables are penalised
with `rho = 10^4` in the objective.

### Changes

**`infNMPC_options.py`** — new field:

| Field | Default | Meaning |
|---|---|---|
| `slack_penalty_weight` | `1e4` | Weight on `m.slack_index` slack vars in the tracking objective. Only applied when `objective='tracking'`; economic objective handles its own penalties via `custom_objective`. Set to 0 to disable. |

**`controllers.py`** — all three tracking objective branches (`FiniteHorizonController` tracking,
`InfiniteHorizonController` Riemann tracking, `InfiniteHorizonController` default tracking) now
append a slack penalty term when `m.slack_index` exists and `slack_penalty_weight > 0`:

```python
stage_cost += slack_penalty_weight * sum(
    getattr(m, sv_name)[idx]
    for sv_name in m.slack_index
    for idx in getattr(m, sv_name).index_set()
)
```

All slack var indices are included (no t=0 exclusion needed — with `clean_model='delete'` and
RADAU, t=0 entries are absent from discretised Vars). For the `InfiniteHorizonController`,
`m.finite_block.slack_index` is used since slack vars live on the finite block.

**`examples/lyap_flag_examples/enmpc_ternary_distillation/run.py`** — weights updated to match
AMPL `double_col_dyn.run` exactly:
- `stage_cost_weights=[10, 10, 10, 1, 1, 1, 1, 1, 1, 1, 1]` — AMPL `state_w=10` for CVs
  (xD1A, xD2B, xC are tray compositions, same category as x1/x2/M state vars in AMPL),
  `cont_w=1` for MVs (VB1, LT1, D1, B1, VB2, LT2, D2, B2)
- `slack_penalty_weight=1e4` — AMPL `rho = 10^4`
- `terminal_soft_weight=1e4` — AMPL `rho * termeps`
- `nfe_finite=3` — shorter horizon than AMPL's 25 FEs (for faster testing)
- `num_horizons=100` — more MPC iterations than AMPL's K=25 (for more results)

---

## `objective` / `tracking_setpoint` Options Refactor (2026-04-21)

**Branch:** `main`

Replaced the boolean `custom_objective` option with two new string options that together
cover the same cases and add a new one.

---

### Motivation

`custom_objective=True/False` conflated two independent decisions:
1. **What is the MPC stage cost?** (economic function vs. quadratic tracking)
2. **Where do the tracking setpoints come from?** (model file vs. economic optimum)

The new options separate these concerns and expose the previously impossible combination:
track a quadratic cost toward the *economically optimal* steady state (useful when the
economic objective defines a natural operating point but the controller should use a
smooth quadratic stage cost rather than the raw economic function).

---

### New options (`infNMPC_options.py`)

| Field | Default | Meaning |
|---|---|---|
| `objective` | `'economic'` | `'economic'` = use model's `custom_objective()` as stage cost; `'tracking'` = quadratic tracking |
| `tracking_setpoint` | `'model'` | Only when `objective='tracking'`: `'model'` = use `m.setpoints`; `'economic'` = solve SS with `custom_objective()`, use result as tracking target |

`custom_objective: bool` is removed entirely.

---

### Mapping from old to new

| Old | New |
|---|---|
| `custom_objective=True` | `objective='economic'` |
| `custom_objective=False` | `objective='tracking'` (default `tracking_setpoint='model'`) |
| *(impossible before)* | `objective='tracking', tracking_setpoint='economic'` |

---

### Code changes

**`infNMPC_options.py`** — replaced `custom_objective` field; added `objective` and
`tracking_setpoint` with full docstrings.

**`make_model.py`**:
- `_solve_steady_state_model`: condition simplified to `if target is None:` — always uses
  economic SS objective when no target is supplied (target is None for both
  `objective='economic'` and `tracking_setpoint='economic'`).
- `_make_infinite_horizon_model` / `_make_finite_horizon_model`: `m_ss_target` is `None`
  when `objective='economic'` or `tracking_setpoint='economic'`; `ScalarData(setpoints)`
  otherwise. `ss_obj_value` is only stored on blocks when `objective='economic'` (the only
  case where it is referenced in the controller objective).
- `_infinite_block_gen._add_dphidt`: phi ODE branches on `options.objective == 'economic'`.

**`controllers.py`**: both `InfiniteHorizonController` and `FiniteHorizonController`
objective branches gate on `options.objective == 'economic'`.

**All run scripts** (9 examples + `live_plot_examples/enmpc_cstr.py` + `tests/run_ternary_lyap.py`)
updated.

---

### `enmpc_ternary_distillation` run.py updated to new combination

The lyap_flag_examples ternary distillation run.py was updated to use
`objective='tracking', tracking_setpoint='economic'` — the new capability.
This finds the economically optimal SS operating point, then runs a finite-horizon
quadratic tracking controller toward it. `lyap_flag=False` in this configuration.

---

## Plant Squareness & Warm-Start (2026-04-15) — RESOLVED

**Branch:** `shifting-behavior`

### Problem: plant solve fails with infeasibility

The ternary distillation plant solve was failing with IPOPT infeasibility on every iteration.
All root causes were identified and fixed.

---

### Root cause 1: Dual infeasibility in controller (FIXED)

`lyap_soft_weight = 1e4` combined with `stage_cost_weights = [1e4, 1e4, 1e4, ...]` produced
Lagrange multipliers of O(1e8) for the Lyapunov constraint, causing `inf_du` to blow up to 1e16
while the point was primal-feasible.

**Fix:** `lyap_soft_weight = 1.0` in `run.py`.

---

### Root cause 2: `load_data` skips fixed variables — MVs and slacks never loaded (FIXED)

`DynamicModelInterface.load_data` skips fixed variables. MVs are fixed in `Plant.__init__` and
slacks are fixed after each MPC iteration. After every iteration, the MV unfix/load/refix loop
in `run_MPC.py` was added so that the plant receives the controller's optimal MV values before
each plant solve. Slacks are similarly unfixed before the warm-start load, then refixed after.

**Fix:** In `run_MPC.py`, before `plant.interface.load_data(plant_warmstart)`:
  - Unfix all MVs (`_mv.unfix()`)
  - Unfix all slacks (`_sv.unfix()`)
  - Call `load_data`
  - Re-fix MVs and slacks

---

### Root cause 3: `lyap_slack` DOF in plant (FIXED)

When `lyap_flag=True`, `Plant.__init__` inherited the lyap options and built `lyap_slack` as a
free `NonNegativeReals` Var with a constant plant objective. This was a spurious DOF.

**Fix:** `plant_options = options.copy(... lyap_flag=False)` in `plant.py`.

---

### Root cause 4: equations_write creates constraints at t=0 — a non-collocation point (FIXED)

This was the deepest root cause, responsible for plant failures at iterations 1+ even after
MVs/slacks were correctly propagated.

**Background:** In LAGRANGE-RADAU with ncp=3, collocation points are at
`{0.155051, 0.644949, 1.0}`. `t=0` is NOT a collocation point. `equations_write` runs AFTER
discretization (necessary because Pyomo's DAE dictionary expansion only works once `m.time` is
populated). As a result, every `Constraint(m.time, rule=...)` in `equations_write` creates an
entry at `t=0`, including algebraic-variable definitions such as:

```python
m.xD2B_def = Constraint(m.time, rule=lambda m, t: m.xD2B[t] == m.x2[41, 2, t])
```

The `xD2B_def[0]` entry actively constrained `x2[41, 2, 0]` (a differential state var) to equal
`xD2B[0]` (an algebraic CV fixed at the initialization value). At iteration 0, the initialization
is consistent so this is harmless. At iteration 1+, the plant IC (state vars at t=0) is updated
to the plant endpoint, but `xD2B[0]` remains stale (fixed at the IC-solve value, never updated
because `load_data` skips fixed vars and the plant has no `shift_values_by_time`). IPOPT then
forced `x2[41, 2, 0]` = the stale `xD2B[0]` value, giving the ODE a **wrong initial condition**
and causing infeasibility at the horizon endpoint.

**Diagnosis method:** Wrote pprints of both controller post-solve and plant pre-solve to
`research_files/shifting_behavior/iter00{0,1}_*.txt`. Comparing `xD2B[0]` between controller
(0.9142, updated by `shift_values_by_time`) and plant (0.8999, stale) revealed the discrepancy.

**Fix — three-part:**

1. **`make_model.py` — `_finite_block_gen`**: After `_model.equations_write(m)`, delete all
   constraint entries at `t=0`:
   ```python
   _t0_del = m.time.first()
   for _con in m.component_objects(pyo.Constraint):
       _keys_at_t0 = [_idx for _idx in list(_con.keys())
                      if (_idx[-1] if isinstance(_idx, tuple) else _idx) == _t0_del]
       for _idx in _keys_at_t0:
           del _con[_idx]
   ```

2. **`make_model.py` — `_fix_initial_conditions_at_t0`**: Simplified entirely. The old function
   fixed algebraic CVs at t=0 and left the state vars they determined (via equality constraints)
   free. With no constraints at t=0, this logic is dead. New version simply fixes **all
   differential state vars at t=0** and nothing else.

3. **`run_MPC.py` — IC loading**: Since all state vars at t=0 are now fixed (no algebraic
   bypass), `load_data` would skip them all. Both the plant and controller IC loads now use an
   explicit unfix → load → refix sequence:
   ```python
   # For plant:
   plant.interface.load_data(tf_data)               # warm-start t>0
   for _sv in plant._model.state_vars: _sv[t0].unfix()
   plant.interface.load_data(tf_data, time_points=t0_plant)
   for _sv in plant._model.state_vars: _sv[t0].fix()

   # For controller: same pattern after shift_values_by_time
   ```

**Result:** All 5 MPC iterations in the ternary distillation lyap example solve optimally.
Plant constraint violations at numerical precision (~1e-12 to 1e-15) only.

---

### Remaining cleanup items (not done)

- `run_MPC.py` line 140–142: unconditional `pprint` to `model_output.txt` at `i==0` (was `i==10`, moved by prior work); not gated by `debug_flag`.
- `plant.py` lines 74–75: unconditional `pprint` to `plant_model.txt` in `Plant.__init__`.
- `run.py`: `num_horizons=5` was set for testing; should be restored to 100 (or intended value) for production runs.
- The `equations_write` fix is only applied in `_finite_block_gen` (RADAU finite block). The infinite block (`_infinite_block_gen`) uses LAGRANGE-LEGENDRE, where BOTH endpoints (`t=0` and `t=t_last`) are non-collocation points. The same spurious constraint issue could exist there. Not yet investigated.

---

## Lyapunov Soft Constraint & Warm-Start Diagnostics (2026-04-13)

**Branch:** `shifting-behavior`

---

### Problem 1: `_ipopt_warm_solver` settings caused convergence to local infeasibility

The initial warm solver implementation used three extra IPOPT options that backfired:

- `warm_start_init_point = 'yes'` — tells IPOPT to use current primal **and dual** (multiplier)
  values as the initial iterate.  Pyomo's `ipopt_v2` interface does not pass dual variable values
  to IPOPT without explicit `Suffix` setup, so IPOPT received zero multipliers.  Starting with
  zero duals and a small barrier produced a badly conditioned initial search direction.
- `mu_init = 1e-3` — fine when the warm-started point is nearly feasible, but the Lyapunov
  constraint violation at the warm-started point grows with each iteration (from `1.27e-7` →
  `5.00e-1` → `1.91e+01` across the first three MPC steps).  With a barrier that tight, IPOPT
  couldn't take large enough steps to recover from `inf_pr = 19.1` and converged to a point of
  local infeasibility after 216 iterations.

**Fix:** removed `warm_start_init_point` and `mu_init` from `_ipopt_warm_solver`.  The primal
warm start is already fully effective — `shift_values_by_time` places the shifted solution
directly in the Pyomo model and IPOPT reads those values as its initial iterate automatically.
`bound_push = 1e-8` / `bound_frac = 1e-8` are retained: they prevent IPOPT from projecting
bound-constrained variables away from the good warm-started positions at the start of each solve.

---

### Problem 2: Lyapunov constraint violation at the warm-started point

Even with the corrected warm solver the failure persisted, because making the warm solver
identical to the cold solver also failed at the same iteration.

**Diagnosis — the doubled endpoint:**

After `shift_values_by_time(sampling_time)` the last time point (e.g. `t = 3` in a `[0, 3]`
horizon) has no `x_old[t + h]` to shift from, so Pyomo keeps the old endpoint value.  The
Riemann-sum `phi_track` then double-counts that endpoint:

```
phi_track_warm ≈ cost_opt(t=2) + cost_opt(t=3) + cost_opt(t=3)   [t=3 repeated]
V_prev          = cost_opt(t=1) + cost_opt(t=2) + cost_opt(t=3)
```

When the optimal trajectory does not have monotonically decreasing stage costs along the
horizon (which the economic + Lyapunov objective does not guarantee), `cost_opt(t=3) >
cost_opt(t=1)` and `phi_track_warm > V_prev`, violating the Lyapunov constraint by the
full amount at the warm-started point.  For `iter 1 → iter 2`, the violation was measured
at exactly `+19.1238`, matching IPOPT's reported `inf_pr = 1.91e+01`.

Running without `lyap_flag` confirmed the physics warm start is healthy (initial `inf_pr`
stays `< 1.2` for all tested iterations).  The Lyapunov constraint is the sole source of
the large infeasibility.

**Root cause:** with a short finite horizon (`nfe_finite = 3`, `sampling_time = 1`), the
controller lacks enough foresight to guarantee Lyapunov decrease while also satisfying the
model physics.  The hard inequality becomes genuinely locally infeasible, and no warm start
tuning can escape a locally infeasible region.

---

### Fix: `lyap_constraint_type` option

Added three-way control over how the Lyapunov decrease condition is enforced, matching the
pattern already used by `terminal_constraint_type`.

**`infNMPC_options.py`** — two new fields:

| Field | Default | Meaning |
|---|---|---|
| `lyap_constraint_type` | `'hard'` | `'hard'` = strict inequality (existing behavior, backward compatible); `'soft'` = L1-relaxed with a non-negative slack; `'none'` = infrastructure built, no constraint |
| `lyap_soft_weight` | `1e4` | Penalty weight on `lyap_slack` when `'soft'` |

**`make_model.py`** — both Lyapunov constraint sites (`_finite_block_gen` for finite-only
controllers, `_make_infinite_horizon_model` for infinite-horizon) now branch on
`lyap_constraint_type`:
- `'hard'`: original `pyo.Constraint` unchanged.
- `'soft'`: adds `m.lyap_slack` (`NonNegativeReals`, initialize=0) and writes
  `phi_track - V_prev <= -lyap_delta * first_stage_cost_prev + lyap_slack`.
- `'none'`: phi_track / V_prev / first_stage_cost_prev infrastructure is still built
  (so `run_MPC.py`'s parameter updates work), but no constraint or variable is added.

**`controllers.py`** — `_use_soft_lyap` flag added to both `InfiniteHorizonController` and
`FiniteHorizonController`.  All five objective rules (`custom`, `riemann`, `default` for
infinite; `custom`, `default` for finite) append
`+ m.lyap_slack * options.lyap_soft_weight` when the flag is set.

**`examples/lyap_flag_examples/enmpc_ternary_distillation/run.py`** — switched to
`lyap_constraint_type='soft'`, `lyap_soft_weight=1e4`.

---

### Also done this session

- **`_ipopt_warm_solver`** added to `make_model.py`; `Controller` now uses the cold solver
  for the initial solve and the warm solver for all subsequent MPC iterations via the
  `_initialized` flag.
- **Shift observation files** added to `run_MPC.py` (gated on `debug_flag and i == 0`):
  `research_files/shifting_behavior/iter000_post_solve.txt` and `iter000_post_shift.txt`.
  `research_files/` is already in `.gitignore`.
- **Typo fix:** `optimisation` → `optimization` in docstrings in `controllers.py` and
  `run_MPC.py`.

---

## Warm-Start Solver for MPC Updates (2026-04-13)

**Branch:** `shifting-behavior`

### Motivation

The MPC loop in `run_MPC.py` already performs a proper shift at each iteration
(lines 245-246):
```python
controller.interface.shift_values_by_time(options.sampling_time)
controller.interface.load_data(tf_data, time_points=t0_controller)
```
This places the shifted previous solution in the model as a natural warm start,
but the controller was using the same cold-start IPOPT settings for every solve.
Warm starts benefit from a smaller initial barrier parameter and tighter bound
handling so IPOPT stays close to the shifted solution.

### Changes

**`infeNMPC/make_model.py`**

- `_ipopt_solver()` — retained as the **cold-start** solver.  Comment placeholders
  for `bound_push` / `mu_init` removed; settings are now explicit and documented
  in the new function.
- `_ipopt_warm_solver()` — **new** function for warm-started MPC updates.  Extra
  settings vs. cold solver:
  - `warm_start_init_point = 'yes'` — IPOPT accepts primal/dual values already
    in the model rather than re-projecting from scratch.
  - `mu_init = 1e-3` — start barrier near the solution rather than the default
    0.1.
  - `bound_push = 1e-8`, `bound_frac = 1e-8` — don't push the starting iterate
    away from bounds; appropriate when it's already feasible.

**`infeNMPC/controllers.py`**

- `Controller.__init__` — instantiates both `_solver` (cold) and `_warm_solver`.
  Adds `_initialized = False` flag.
- `Controller.solve()` — selects solver based on `_initialized`:
  - First call (from `__init__` via subclass): uses cold solver.
  - All subsequent calls (MPC loop after shift): uses warm solver.
  Sets `_initialized = True` at end of first successful solve.

---

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

## ~~Fix initial condition solve — always solve the iv model~~ (Done)

**Done as of binary-distillation branch.** The `if all(var in CV_index...)` shortcut no longer
exists in `make_model.py`. Both `_make_finite_horizon_model` and `_make_infinite_horizon_model`
call `_build_ic_data(options)` unconditionally, which routes to `_build_ic_data_steady_state` →
`_solve_steady_state_model`. The returned `ScalarData` contains consistent values for all
variables (not just CVs), so every state variable at `t=0` is correctly initialized.

**Original description of the bug and fix (preserved for history):**

Both orchestrators formerly contained:
```python
if all(var in m_iv.CV_index for var in initial_value_vars):
    initial_data = ScalarData({...})   # SKIP — only CV values, no state vars
else:
    initial_data = _solve_steady_state_model(...)
```
When `initial_values` only contained CVs (every current example), the iv solve was skipped,
leaving the 240+ non-CV state variables at their potentially inconsistent `initialize=` guesses.
Fix: remove the branch; always call `_solve_steady_state_model` via `_build_ic_data`.

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
