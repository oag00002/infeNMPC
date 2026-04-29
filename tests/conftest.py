"""
Shared pytest fixtures for the infeNMPC test suite.

All integration tests use the eNMPC CSTR model (examples/lyap_flag_examples/
enmpc_cstr/small_cstr_model.py) — it is the smallest, fastest-converging model
in the repository and is the only model confirmed working unconditionally.

Fixtures are parameterised so that different option combinations can be tested
without duplicating the model-import boilerplate.
"""
import sys
from pathlib import Path

import pytest

# ---------------------------------------------------------------------------
# Make the CSTR example model importable from tests
# ---------------------------------------------------------------------------
_CSTR_DIR = (
    Path(__file__).parent.parent
    / "examples"
    / "lyap_flag_examples"
    / "enmpc_cstr"
)
if str(_CSTR_DIR) not in sys.path:
    sys.path.insert(0, str(_CSTR_DIR))


@pytest.fixture(scope="session")
def cstr_module():
    """Return the CSTR model module (imported once per session)."""
    import small_cstr_model  # noqa: PLC0415
    return small_cstr_model


@pytest.fixture
def base_options(cstr_module):
    """
    Minimal Options for the CSTR that completes in one or two IPOPT solves.

    num_horizons=1 means a single controller solve + single plant integrate.
    nfe=1 / ncp=1 gives the coarsest legal discretization — fast, but still
    exercises the full model-construction path.
    """
    from infeNMPC import Options
    return Options.for_model_module(
        cstr_module,
        num_horizons=1,
        sampling_time=1.0,
        nfe_finite=1,
        ncp_finite=1,
        nfe_infinite=1,
        ncp_infinite=1,
        infinite_horizon=True,
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


@pytest.fixture
def lyap_options(cstr_module):
    """
    Options for the CSTR with Lyapunov stability constraint enabled.

    Uses the same discretization as the canonical run.py but with only
    3 MPC iterations so tests finish quickly.
    """
    from infeNMPC import Options
    return Options.for_model_module(
        cstr_module,
        num_horizons=3,
        sampling_time=1.0,
        nfe_finite=2,
        ncp_finite=1,
        nfe_infinite=5,
        ncp_infinite=1,
        infinite_horizon=True,
        terminal_constraint_type="none",
        objective="economic",
        stage_cost_weights=[1, 1, 1 / 600],
        beta=1.0,
        lyap_flag=True,
        lyap_delta=0.01,
        lyap_constraint_type="hard",
        save_data=False,
        save_figure=False,
        tee_flag=False,
        debug_flag=False,
    )
