"""
infeNMPC command-line entry point.

Usage
-----
Run an example by pointing at its model file:

    python -m infeNMPC --model examples/enmpc_cstr/model.py
    python -m infeNMPC --model examples/pendulum/model.py

Or run an example's ``run.py`` directly (preferred for non-default settings):

    python examples/enmpc_cstr/run.py
    python examples/pendulum/run.py
"""
import argparse
import importlib.util
import sys
from pathlib import Path

from .infNMPC_options import Options
from .run_MPC import mpc_loop


def _load_module_from_path(path: str):
    """Load a Python module from an arbitrary file path."""
    p = Path(path).resolve()
    spec = importlib.util.spec_from_file_location(p.stem, p)
    module = importlib.util.module_from_spec(spec)
    # Make the module's directory importable (so relative imports inside work)
    sys.path.insert(0, str(p.parent))
    spec.loader.exec_module(module)
    return module


def main():
    parser = argparse.ArgumentParser(
        description='Run an infeNMPC closed-loop simulation.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        '--model', '-m',
        required=True,
        metavar='MODEL_FILE',
        help=(
            'Path to a model.py file implementing variables_initialize and '
            'equations_write (and optionally custom_objective and '
            'default_options).  Example: examples/enmpc_cstr/model.py'
        ),
    )
    args = parser.parse_args()

    model_module = _load_module_from_path(args.model)
    options = Options.for_model_module(model_module)
    mpc_loop(options)


if __name__ == '__main__':
    main()
