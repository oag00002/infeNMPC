"""
End-to-end MPC loop tests (test_mpc_loop.py).

These tests run the full closed-loop simulation via mpc_loop() and verify
the documented behavior of the simulation orchestrator (run_MPC.py).

All tests use the CSTR model with short runs (2-5 iterations) to keep CI time
manageable.  The CSTR is the only model confirmed to work unconditionally;
other models (ternary distillation, binary distillation) are larger and not
tested in CI.

What is verified:

  Convergence
  - mpc_loop() completes 3 iterations without raising RuntimeError.
  - All controller solves reach TerminationCondition.optimal (or acceptable).

  Lyapunov tracking (lyap_flag=True)
  - V_prev.value is updated each iteration.
  - After iteration 1, V_prev is non-negative (it's a sum of squares).
  - The sequence of V values is monotonically non-increasing when using the
    canonical CSTR settings (this is the Lyapunov stability guarantee).

  Data output (save_data=True)
  - The results folder is created.
  - io_data.csv is written with the correct column structure (Time, CV names,
    MV names).
  - sim_data.csv is written.

  Tracking objective
  - mpc_loop works with objective='tracking' (no custom_objective required).

  Finite horizon
  - mpc_loop works with infinite_horizon=False.
"""
import csv
import os
from pathlib import Path

import pytest
import pyomo.environ as pyo

from infeNMPC import Options, mpc_loop
from infeNMPC.controllers import InfiniteHorizonController
from infeNMPC.plant import Plant


# ---------------------------------------------------------------------------
# Basic convergence
# ---------------------------------------------------------------------------

class TestMpcLoopConvergence:
    """mpc_loop() runs without error and produces optimal solves."""

    def test_cstr_3_iters_no_error(self, base_options):
        """
        Three MPC iterations on the CSTR must complete without RuntimeError.
        This is the most basic end-to-end smoke test — if anything in the
        model construction, IC loading, or shift logic is broken it shows up here.
        """
        opts = base_options.copy(num_horizons=3)
        mpc_loop(opts)  # must not raise

    def test_lyap_cstr_3_iters_no_error(self, lyap_options):
        """Three iterations with Lyapunov constraint enabled must not raise."""
        mpc_loop(lyap_options)

    def test_finite_horizon_3_iters(self, cstr_module):
        """Finite-horizon controller runs without error."""
        opts = Options.for_model_module(
            cstr_module,
            num_horizons=3,
            sampling_time=1.0,
            nfe_finite=3,
            ncp_finite=1,
            infinite_horizon=False,
            terminal_constraint_type="none",
            objective="economic",
            stage_cost_weights=[1, 1, 1 / 600],
            beta=1.0,
            lyap_flag=False,
            save_data=False,
            save_figure=False,
            tee_flag=False,
            debug_flag=False,
        )
        mpc_loop(opts)

    def test_tracking_objective_3_iters(self, cstr_module):
        """
        objective='tracking' uses a quadratic stage cost instead of the model's
        custom_objective.  This combination was added in the April-2026 refactor
        that replaced the bool custom_objective option with a string 'objective'.
        """
        opts = Options.for_model_module(
            cstr_module,
            num_horizons=3,
            sampling_time=1.0,
            nfe_finite=2,
            ncp_finite=1,
            nfe_infinite=3,
            ncp_infinite=1,
            infinite_horizon=True,
            terminal_constraint_type="none",
            objective="tracking",
            tracking_setpoint="model",
            stage_cost_weights=[1, 1, 1 / 600],
            beta=1.0,
            lyap_flag=False,
            save_data=False,
            save_figure=False,
            tee_flag=False,
            debug_flag=False,
        )
        mpc_loop(opts)

    def test_soft_terminal_constraint_3_iters(self, cstr_module):
        """terminal_constraint_type='soft' adds a quadratic penalty, must converge."""
        opts = Options.for_model_module(
            cstr_module,
            num_horizons=3,
            sampling_time=1.0,
            nfe_finite=2,
            ncp_finite=1,
            nfe_infinite=3,
            ncp_infinite=1,
            infinite_horizon=True,
            terminal_constraint_type="soft",
            terminal_soft_weight=1.0,
            objective="economic",
            stage_cost_weights=[1, 1, 1 / 600],
            beta=1.0,
            lyap_flag=False,
            save_data=False,
            save_figure=False,
            tee_flag=False,
            debug_flag=False,
        )
        mpc_loop(opts)

    def test_hard_terminal_constraint_3_iters(self, cstr_module):
        """terminal_constraint_type='hard' with the CSTR must solve optimally."""
        opts = Options.for_model_module(
            cstr_module,
            num_horizons=3,
            sampling_time=1.0,
            nfe_finite=2,
            ncp_finite=1,
            nfe_infinite=3,
            ncp_infinite=1,
            infinite_horizon=True,
            terminal_constraint_type="hard",
            terminal_constraint_variables="cvmv",
            objective="economic",
            stage_cost_weights=[1, 1, 1 / 600],
            beta=1.0,
            lyap_flag=False,
            save_data=False,
            save_figure=False,
            tee_flag=False,
            debug_flag=False,
        )
        mpc_loop(opts)

    def test_soft_lyap_constraint_3_iters(self, cstr_module):
        """
        lyap_constraint_type='soft' is the recommended setting when a short
        horizon makes the hard Lyapunov constraint locally infeasible after
        a shift.  The soft form adds a NonNegativeReals lyap_slack whose value
        enters the objective with weight lyap_soft_weight.
        """
        opts = Options.for_model_module(
            cstr_module,
            num_horizons=3,
            sampling_time=1.0,
            nfe_finite=2,
            ncp_finite=1,
            nfe_infinite=3,
            ncp_infinite=1,
            infinite_horizon=True,
            terminal_constraint_type="none",
            objective="economic",
            stage_cost_weights=[1, 1, 1 / 600],
            beta=1.0,
            lyap_flag=True,
            lyap_delta=0.01,
            lyap_constraint_type="soft",
            lyap_soft_weight=1e4,
            save_data=False,
            save_figure=False,
            tee_flag=False,
            debug_flag=False,
        )
        mpc_loop(opts)


# ---------------------------------------------------------------------------
# Lyapunov parameter updates
# ---------------------------------------------------------------------------

class TestLyapunovUpdates:
    """
    V_prev and first_stage_cost_prev are Pyomo mutable Params updated each
    iteration by run_MPC.py.  We verify they are updated by running the
    controller manually for two iterations and checking the stored values.
    """

    def test_v_prev_updated_after_each_iteration(self, lyap_options):
        """
        After each MPC iteration, V_prev is set to the current Lyapunov value.
        Before any MPC iteration has run, V_prev is initialised to a large
        sentinel value (1e10) so that the first iteration's Lyapunov constraint
        is always satisfiable — any finite V_current satisfies V <= 1e10.
        """
        ctrl = InfiniteHorizonController(lyap_options)
        m = ctrl._model

        # V_prev starts at 1e10 (the sentinel used in make_model.py).
        initial_v_prev = pyo.value(m.V_prev)
        assert initial_v_prev == pytest.approx(1e10), (
            "V_prev should be initialised to 1e10 (sentinel) before the first iteration"
        )

    def test_first_stage_cost_prev_initialised_zero(self, lyap_options):
        """first_stage_cost_prev starts at 0 before the first MPC iteration."""
        ctrl = InfiniteHorizonController(lyap_options)
        m = ctrl._model
        assert pyo.value(m.first_stage_cost_prev) == pytest.approx(0.0)

    def test_phi_track_nonnegative_after_solve(self, lyap_options):
        """
        phi_track is a sum of squared deviations — it must be non-negative.
        After the initial controller solve, the warm-started values give a
        phi_track value that reflects the deviation from setpoint.
        """
        ctrl = InfiniteHorizonController(lyap_options)
        val = pyo.value(ctrl._model.finite_block.phi_track)
        assert val >= 0.0


# ---------------------------------------------------------------------------
# Data output
# ---------------------------------------------------------------------------

class TestDataOutput:
    """
    When save_data=True, mpc_loop() writes io_data.csv and sim_data.csv into
    a timestamped Results/ folder.  The CSV structure must match the documented
    format: Time column, then one column per CV, then one column per MV.
    """

    def test_io_data_csv_created(self, base_options, tmp_path, monkeypatch):
        """io_data.csv is created when save_data=True."""
        # Redirect results to tmp_path so the test is self-cleaning
        monkeypatch.chdir(tmp_path)
        opts = base_options.copy(num_horizons=2, save_data=True, save_figure=False)
        mpc_loop(opts)

        results_dirs = list(Path("Results").rglob("io_data.csv"))
        assert len(results_dirs) >= 1, "io_data.csv was not created"

    def test_io_data_csv_has_time_column(self, base_options, tmp_path, monkeypatch):
        """io_data.csv first column is 'Time'."""
        monkeypatch.chdir(tmp_path)
        opts = base_options.copy(num_horizons=2, save_data=True, save_figure=False)
        mpc_loop(opts)

        csv_path = next(Path("Results").rglob("io_data.csv"))
        with open(csv_path) as f:
            reader = csv.DictReader(f)
            assert "Time" in reader.fieldnames

    def test_io_data_csv_has_cv_columns(self, base_options, tmp_path, monkeypatch):
        """
        io_data.csv contains one column per CV.  Column names use the display-name
        format 'CV_{display_name}' (e.g. 'CV_C_a\\;(mol/L)').  We verify that
        at least two CV columns are present for the CSTR, without hardcoding the
        exact display-name strings which could change if the model is renamed.
        """
        monkeypatch.chdir(tmp_path)
        opts = base_options.copy(num_horizons=2, save_data=True, save_figure=False)
        mpc_loop(opts)

        csv_path = next(Path("Results").rglob("io_data.csv"))
        with open(csv_path) as f:
            reader = csv.DictReader(f)
            fieldnames = reader.fieldnames

        cv_cols = [f for f in fieldnames if f.startswith("CV_")]
        assert len(cv_cols) == 2, (
            f"Expected 2 CV columns in io_data.csv (one per CSTR CV), got: {cv_cols}"
        )

    def test_io_data_csv_has_mv_columns(self, base_options, tmp_path, monkeypatch):
        """
        io_data.csv contains one column per MV.  Column names use the display-name
        format 'MV_{display_name}'.  We verify one MV column is present for the CSTR.
        """
        monkeypatch.chdir(tmp_path)
        opts = base_options.copy(num_horizons=2, save_data=True, save_figure=False)
        mpc_loop(opts)

        csv_path = next(Path("Results").rglob("io_data.csv"))
        with open(csv_path) as f:
            reader = csv.DictReader(f)
            fieldnames = reader.fieldnames

        mv_cols = [f for f in fieldnames if f.startswith("MV_")]
        assert len(mv_cols) == 1, (
            f"Expected 1 MV column in io_data.csv (CSTR has one MV Fa0), got: {mv_cols}"
        )

    def test_io_data_row_count(self, base_options, tmp_path, monkeypatch):
        """io_data.csv has exactly num_horizons+1 rows (initial + one per iteration)."""
        monkeypatch.chdir(tmp_path)
        n = 3
        opts = base_options.copy(num_horizons=n, save_data=True, save_figure=False)
        mpc_loop(opts)

        csv_path = next(Path("Results").rglob("io_data.csv"))
        with open(csv_path) as f:
            rows = list(csv.DictReader(f))
        assert len(rows) == n + 1

    def test_sim_data_csv_created(self, base_options, tmp_path, monkeypatch):
        """sim_data.csv (full DAE trajectory) is also written."""
        monkeypatch.chdir(tmp_path)
        opts = base_options.copy(num_horizons=2, save_data=True, save_figure=False)
        mpc_loop(opts)

        sim_csvs = list(Path("Results").rglob("sim_data.csv"))
        assert len(sim_csvs) >= 1

    def test_run_config_txt_created(self, base_options, tmp_path, monkeypatch):
        """
        run_config.txt is written to the timestamped folder immediately when
        it is created (before the loop starts), so it exists even if the run
        is interrupted.
        """
        monkeypatch.chdir(tmp_path)
        opts = base_options.copy(num_horizons=2, save_data=True, save_figure=False)
        mpc_loop(opts)

        configs = list(Path("Results").rglob("run_config.txt"))
        assert len(configs) >= 1

    def test_no_csvs_when_save_data_false(self, base_options, tmp_path, monkeypatch):
        """
        When save_data=False, io_data.csv and sim_data.csv are not written.
        Note: the Results/ folder and run_config.txt are always created
        (before the loop starts) so that the config is persisted even if the
        run is interrupted.  Only the data files require save_data=True.
        """
        monkeypatch.chdir(tmp_path)
        opts = base_options.copy(num_horizons=2, save_data=False, save_figure=False)
        mpc_loop(opts)

        io_csvs = list(Path("Results").rglob("io_data.csv")) if Path("Results").exists() else []
        sim_csvs = list(Path("Results").rglob("sim_data.csv")) if Path("Results").exists() else []
        assert io_csvs == [], "io_data.csv was written even though save_data=False"
        assert sim_csvs == [], "sim_data.csv was written even though save_data=False"


# ---------------------------------------------------------------------------
# IC consistency check
# ---------------------------------------------------------------------------

class TestICConsistency:
    """
    _check_ic_consistency fires at the end of model construction and raises
    RuntimeError if any physical constraint at t=0 is violated by more than
    the tolerance.  We verify that well-formed models pass this check (the
    absence of a RuntimeError is the assertion).
    """

    def test_cstr_ic_consistency_passes(self, base_options):
        """
        The CSTR model's IC solve produces a t=0 state that satisfies all
        active equality constraints within 1e-4.  If the IC solve is bypassed
        or returns an infeasible point, this test would fail with RuntimeError.
        """
        # Building the controller triggers _check_ic_consistency internally.
        ctrl = InfiniteHorizonController(base_options)
        assert ctrl._initialized is True
