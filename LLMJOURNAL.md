# Next Tasks

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
