"""
Tests for the Options dataclass (infeNMPC/infNMPC_options.py).

Options is the single configuration object that controls every aspect of an
infeNMPC simulation.  These tests verify:
  - Every field has a documented default value that constructs without error.
  - The for_model_module() classmethod correctly wires up model_module.
  - copy() returns a new instance with overridden fields while leaving the
    original unchanged (immutable-style update pattern).
  - Key fields introduced during refactors (objective, lyap_flag, etc.) are
    present with the correct defaults, guarding against accidental regressions.
"""
import pytest
from infeNMPC import Options


# ---------------------------------------------------------------------------
# Basic instantiation
# ---------------------------------------------------------------------------

class TestOptionsDefaults:
    """Options() with no arguments uses documented defaults."""

    def test_default_construction(self):
        """Options() must not raise — every field has a default."""
        opts = Options()
        assert opts is not None

    def test_infinite_horizon_default_true(self):
        """Infinite-horizon MPC is enabled by default."""
        assert Options().infinite_horizon is True

    def test_nfe_finite_default(self):
        """Default finite-block discretization: 2 finite elements."""
        assert Options().nfe_finite == 2

    def test_ncp_finite_default(self):
        """Default finite-block collocation: 3 points per element."""
        assert Options().ncp_finite == 3

    def test_nfe_infinite_default(self):
        """Default infinite-block discretization: 3 finite elements."""
        assert Options().nfe_infinite == 3

    def test_ncp_infinite_default(self):
        """Default infinite-block collocation: 3 points per element."""
        assert Options().ncp_infinite == 3

    def test_sampling_time_default(self):
        """Default sampling interval is 1.0 time unit."""
        assert Options().sampling_time == 1.0

    def test_num_horizons_default(self):
        """Default run length is 100 MPC iterations."""
        assert Options().num_horizons == 100

    def test_objective_default_economic(self):
        """
        objective='economic' is the default, meaning the model's custom_objective()
        is used as the stage cost.  This changed from a bool (custom_objective=True)
        in the April-2026 refactor; the string form is checked here so a regression
        to the old bool API is immediately caught.
        """
        assert Options().objective == "economic"

    def test_tracking_setpoint_default(self):
        """tracking_setpoint='model' uses m.setpoints from the model file."""
        assert Options().tracking_setpoint == "model"

    def test_terminal_constraint_type_default(self):
        """
        Terminal constraints default to 'hard' (equality at terminal time).
        The three-way string option replaced the old bool endpoint_constraints
        in the April-2026 overhaul.
        """
        assert Options().terminal_constraint_type == "hard"

    def test_terminal_constraint_variables_default(self):
        """Both CVs and MVs are constrained at the terminal time by default."""
        assert Options().terminal_constraint_variables == "cvmv"

    def test_lyap_flag_default_false(self):
        """Lyapunov infrastructure is disabled by default."""
        assert Options().lyap_flag is False

    def test_lyap_constraint_type_default(self):
        """Lyapunov constraint is hard by default when lyap_flag=True."""
        assert Options().lyap_constraint_type == "hard"

    def test_lyap_delta_default(self):
        """lyap_delta=0.01 requires only a 1% Lyapunov decrease per step."""
        assert Options().lyap_delta == 0.01

    def test_lyap_beta_default(self):
        """lyap_beta=1.2 scales the infinite-block contribution to V."""
        assert Options().lyap_beta == 1.2

    def test_revive_run_default_zero(self):
        """
        revive_run=0 disables the retry-on-failure logic added in April-2026
        to handle 'restoration phase converged to small primal infeasibility'.
        """
        assert Options().revive_run == 0

    def test_save_data_default_true(self):
        """Results CSVs are saved by default."""
        assert Options().save_data is True

    def test_dynamic_initial_conditions_default_false(self):
        """
        dynamic_initial_conditions=False means the IC solve uses a steady-state
        NLP.  The dynamic (non-SS) path must be explicitly opted into.
        """
        assert Options().dynamic_initial_conditions is False

    def test_gamma_default_none(self):
        """
        gamma=None triggers auto-computation from the collocation scheme.
        Passing an explicit float overrides the auto-computation.
        """
        assert Options().gamma is None

    def test_model_module_default_none(self):
        """model_module=None until set by for_model_module()."""
        assert Options().model_module is None


# ---------------------------------------------------------------------------
# for_model_module() classmethod
# ---------------------------------------------------------------------------

class TestForModelModule:
    """for_model_module() wires up the live model object and overrides."""

    def test_sets_model_module(self, cstr_module):
        """model_module is stored on the returned Options."""
        opts = Options.for_model_module(cstr_module)
        assert opts.model_module is cstr_module

    def test_overrides_applied(self, cstr_module):
        """Keyword overrides are applied on top of defaults."""
        opts = Options.for_model_module(cstr_module, num_horizons=7)
        assert opts.num_horizons == 7

    def test_model_name_set_from_module(self, cstr_module):
        """model_name is inferred from the module's __name__."""
        opts = Options.for_model_module(cstr_module)
        assert opts.model_name != "unknown"

    def test_multiple_overrides_applied(self, cstr_module):
        """Multiple keyword overrides are all applied."""
        opts = Options.for_model_module(
            cstr_module, num_horizons=7, sampling_time=2.0, lyap_flag=True
        )
        assert opts.num_horizons == 7
        assert opts.sampling_time == 2.0
        assert opts.lyap_flag is True


# ---------------------------------------------------------------------------
# copy() immutable-style update
# ---------------------------------------------------------------------------

class TestOptionsCopy:
    """copy() returns a new Options with selective field overrides."""

    def test_copy_returns_new_object(self, cstr_module):
        opts = Options.for_model_module(cstr_module)
        opts2 = opts.copy(num_horizons=5)
        assert opts2 is not opts

    def test_copy_overrides_field(self, cstr_module):
        opts = Options.for_model_module(cstr_module, num_horizons=10)
        opts2 = opts.copy(num_horizons=5)
        assert opts2.num_horizons == 5

    def test_copy_does_not_mutate_original(self, cstr_module):
        opts = Options.for_model_module(cstr_module, num_horizons=10)
        opts.copy(num_horizons=5)
        assert opts.num_horizons == 10

    def test_copy_preserves_unspecified_fields(self, cstr_module):
        """Fields not mentioned in copy() retain their original values."""
        opts = Options.for_model_module(cstr_module, sampling_time=2.0, num_horizons=10)
        opts2 = opts.copy(num_horizons=5)
        assert opts2.sampling_time == 2.0
