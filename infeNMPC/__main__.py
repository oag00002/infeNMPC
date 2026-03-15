"""
infeNMPC command-line entry point.

Usage
-----
Run with default settings from ``infNMPC_options.py``:

    python -m infeNMPC

Override the model at runtime (all other settings still come from
``infNMPC_options.py``):

    python -m infeNMPC --model enmpc_cstr
    python -m infeNMPC --model pendulum
    python -m infeNMPC --model nonisothermal_cstr
    python -m infeNMPC --model binary_distillation
    python -m infeNMPC --model enmpc_binary_distillation

Optionally apply the model's built-in default options instead of the
values in ``infNMPC_options.py``:

    python -m infeNMPC --model pendulum --use-model-defaults
"""
import argparse
from .infNMPC_options import _import_settings
from .run_MPC import _mpc_loop
from .model_equations import _load_model


def main():
    parser = argparse.ArgumentParser(
        description='Run an infeNMPC closed-loop simulation.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        '--model', '-m',
        default=None,
        metavar='MODEL_NAME',
        help=(
            'Model to simulate.  Must match a module in infeNMPC/models/. '
            'E.g.: enmpc_cstr, pendulum, nonisothermal_cstr, '
            'binary_distillation, enmpc_binary_distillation. '
            'Defaults to the model_name set in infNMPC_options.py.'
        ),
    )
    parser.add_argument(
        '--use-model-defaults',
        action='store_true',
        default=False,
        help=(
            'If set, apply the selected model\'s default_options() overrides '
            'on top of infNMPC_options.py values.'
        ),
    )
    args = parser.parse_args()

    options = _import_settings()

    if args.model is not None:
        options.model_name = args.model

    if args.use_model_defaults:
        model_module = _load_model(options.model_name)
        if hasattr(model_module, 'default_options'):
            for key, val in model_module.default_options().items():
                setattr(options, key, val)

    _mpc_loop(options)


if __name__ == '__main__':
    main()
