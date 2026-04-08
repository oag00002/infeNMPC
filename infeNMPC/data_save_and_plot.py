import os
import csv
import re
import subprocess
from dataclasses import fields
from datetime import datetime

import pyomo.contrib.mpc as mpc
import matplotlib.pyplot as plt
import numpy as np


# ── Plot script template ───────────────────────────────────────────────────────
# Filled in by _write_plot_script().  All literal {{ / }} are escaped braces in
# the generated Python source; {PLACEHOLDER} tokens are substitution sites.
_PLOT_SCRIPT_TEMPLATE = '''\
#!/usr/bin/env python3
"""
plot_results.py  —  auto-generated plotting script
  Model     : {model_name}
  Run date  : {run_datetime}
  Git branch: {git_branch}

Edit the CONFIGURATION section below, then run:
    python plot_results.py
Figures are written to  figures/  within this directory.
Both PDF (vector, for papers) and PNG (raster preview) are produced.
"""
import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

# ── MATPLOTLIB / LATEX STYLE ──────────────────────────────────────────────────
# Set text.usetex = True if a full LaTeX installation is available on your PATH;
# this enables real LaTeX rendering (required for some fonts / packages).
# With text.usetex = False, matplotlib\'s built-in mathtext handles $...$ labels.
mpl.rcParams.update({{
    "text.usetex"    : False,
    "font.family"    : "serif",
    "font.size"      : 10,
    "axes.titlesize" : 11,
    "axes.labelsize" : 10,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "lines.linewidth": 1.5,
    "figure.dpi"     : 150,    # screen preview
    "savefig.dpi"    : 300,    # saved files
    "savefig.bbox"   : "tight",
}})

# ── CONFIGURATION ─────────────────────────────────────────────────────────────

# Each entry produces one figure:
#   (csv_column, y_axis_label, figure_title)
# csv_column must exactly match a header in sim_data.csv (or io_data.csv).
# Labels / titles are mathtext strings — wrap math in $...$
PLOTS = [
{plots_block}
]

# Horizontal setpoint / reference lines per column: {{csv_column: value}}
SETPOINTS = {{
{setpoints_block}
}}

# Figure dimensions (inches).  3.5" = single journal column; 7.0" = double column.
FIGSIZE = (3.5, 2.8)

# x-axis label
TIME_LABEL = {time_label_repr}

# Source CSV — switch to "io_data.csv" for coarser sampling-interval data
DATA_FILE = "sim_data.csv"

FIGURES_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")

# ── END CONFIGURATION ─────────────────────────────────────────────────────────


def _safe_name(col: str) -> str:
    """Convert a csv column key like \'Ca[*]\' to a filename-safe string."""
    return (col.replace("[", "_").replace("]", "")
               .replace(",", "_").replace("*", "all").strip("_"))


def main() -> None:
    data_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), DATA_FILE)
    df = pd.read_csv(data_path)
    os.makedirs(FIGURES_DIR, exist_ok=True)

    for col, y_label, title in PLOTS:
        if col not in df.columns:
            print(f"[warn] column {{col!r}} not found in {{DATA_FILE}} — skipping")
            continue

        fig, ax = plt.subplots(figsize=FIGSIZE)
        ax.plot(df["Time"], df[col], color="#1f77b4", label="Trajectory")

        if col in SETPOINTS:
            ax.axhline(SETPOINTS[col], color="#d62728", linestyle="--",
                       linewidth=1.0, label="Setpoint")
            ax.legend(framealpha=0.7)

        ax.set_xlabel(TIME_LABEL)
        ax.set_ylabel(y_label)
        ax.set_title(title)
        ax.grid(True, linewidth=0.4, alpha=0.6)
        fig.tight_layout()

        name = _safe_name(col)
        for ext in ("pdf", "png"):
            out = os.path.join(FIGURES_DIR, f"{{name}}.{{ext}}")
            fig.savefig(out)
            print(f"  Saved {{out}}")
        plt.close(fig)

    print(f"\\nDone. Figures written to: {{FIGURES_DIR}}")


if __name__ == "__main__":
    main()
'''


def _write_plot_script(folder_path, options, plant, sim_data):
    """
    Write a customised ``plot_results.py`` into *folder_path*.

    The script is pre-configured with:
    * Active PLOTS entries for every CV and MV (using their display names as
      LaTeX y-axis labels).
    * Commented-out entries for every other column found in sim_data (easy to
      uncomment for state-variable diagnostics).
    * SETPOINTS filled from plant.steady_state_values for each CV.

    Parameters
    ----------
    folder_path : str
        Destination directory (the timestamped run folder).
    options : Options
        Simulation configuration (provides model_name).
    plant : Plant
        Provides CV/MV index lists, display names, and steady-state values.
    sim_data : TimeSeriesData
        Used only to enumerate all available column keys.
    """
    from .tools.indexing_tools import _get_variable_key_for_data

    # ── Gather column names ────────────────────────────────────────────────────
    all_sim_cols = list(sim_data.get_data().keys())

    cv_index = list(plant.CV_index)
    mv_index = list(plant.MV_index)
    cv_display = list(plant.CV_display_names)
    mv_display = list(plant.MV_display_names)

    cv_cols = [_get_variable_key_for_data(plant, v) for v in cv_index]
    mv_cols = [_get_variable_key_for_data(plant, v) for v in mv_index]
    active_cols = set(cv_cols) | set(mv_cols)

    # ── Build PLOTS block ──────────────────────────────────────────────────────
    plot_lines = []

    plot_lines.append("    # ── Controlled Variables ──────────────────────────────────────────────")
    for col, display in zip(cv_cols, cv_display):
        y_label = f"${display}$"
        title   = f"${display}$"
        plot_lines.append(f"    ({col!r}, {y_label!r}, {title!r}),")

    plot_lines.append("    # ── Manipulated Variables ─────────────────────────────────────────────")
    for col, display in zip(mv_cols, mv_display):
        y_label = f"${display}$"
        title   = f"${display}$"
        plot_lines.append(f"    ({col!r}, {y_label!r}, {title!r}),")

    other_cols = [c for c in all_sim_cols if c not in active_cols]
    if other_cols:
        plot_lines.append(
            "    # ── Other state variables — uncomment to plot ────────────────────────"
        )
        for col in other_cols:
            placeholder_label = f"${col}$"
            plot_lines.append(
                f"    # ({col!r}, {placeholder_label!r}, {placeholder_label!r}),"
            )

    plots_block = "\n".join(plot_lines)

    # ── Build SETPOINTS block ──────────────────────────────────────────────────
    sp_lines = []
    ss_vals = getattr(plant, "steady_state_values", {})
    for col, var_name in zip(cv_cols, cv_index):
        val = ss_vals.get(var_name)
        if val is not None:
            sp_lines.append(f"    {col!r}: {val!r},")
    setpoints_block = "\n".join(sp_lines)

    # ── Time label ────────────────────────────────────────────────────────────
    time_display = getattr(plant, "time_display_name", ["Time"])
    time_label_repr = repr(time_display[0]) if time_display else repr("Time")

    # ── Fill template and write ────────────────────────────────────────────────
    now_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    script = _PLOT_SCRIPT_TEMPLATE.format(
        model_name=options.model_name,
        run_datetime=now_str,
        git_branch=_get_git_branch(),
        plots_block=plots_block,
        setpoints_block=setpoints_block,
        time_label_repr=time_label_repr,
    )

    script_path = os.path.join(folder_path, "plot_results.py")
    with open(script_path, "w") as f:
        f.write(script)
    print(f"[Info] Plotting script written to {script_path}")


def _get_git_branch():
    """Return the current git branch name, or 'unknown' if not in a git repo."""
    try:
        result = subprocess.run(
            ["git", "rev-parse", "--abbrev-ref", "HEAD"],
            capture_output=True, text=True, timeout=5,
        )
        return result.stdout.strip() or "unknown"
    except Exception:
        return "unknown"


def _create_run_folder(options, resolved_gamma=None):
    """
    Create and return the results directory for this run.

    The path structure is ``Results/{model_name}/{YYYY-MM-DD}/{HH-MM-SS}/``.
    A ``run_config.txt`` file is written immediately into the new folder
    containing all ``Options`` field values (with ``gamma`` replaced by the
    resolved value when auto-computed) and the current git branch.

    Parameters
    ----------
    options : Options
        Simulation configuration for the run.
    resolved_gamma : float or None
        The gamma value actually used (after auto-computation).  Overrides
        ``options.gamma`` in the config file when not None.

    Returns
    -------
    str
        Path to the newly created results directory.
    """
    now = datetime.now()
    date_str = now.strftime("%Y-%m-%d")
    time_str = now.strftime("%H-%M-%S")

    model_name = re.sub(r'[^a-zA-Z0-9_\-]', '_', options.model_name or 'unknown')
    folder_path = os.path.join("Results", model_name, date_str, time_str)
    os.makedirs(folder_path, exist_ok=True)

    # Write run_config.txt
    config_path = os.path.join(folder_path, "run_config.txt")
    git_branch = _get_git_branch()
    with open(config_path, "w") as f:
        f.write(f"git_branch: {git_branch}\n\n")
        f.write("[options]\n")
        for fld in fields(options):
            if fld.name == "model_module":
                continue  # not serialisable
            value = getattr(options, fld.name)
            if fld.name == "gamma" and resolved_gamma is not None:
                value = resolved_gamma
            f.write(f"{fld.name}: {value}\n")

    return folder_path


def _get_results_folder(options):
    """
    Legacy helper: return a results folder path without a timestamp.

    New code should use ``_create_run_folder`` instead.  This function is
    retained for ``_save_epsilon`` and ``_save_lyap_csv`` which may be called
    outside the main MPC loop and therefore have no run-start timestamp.

    The path structure is ``Results/{model_name}/``.
    """
    model_name = re.sub(r'[^a-zA-Z0-9_\-]', '_', options.model_name or 'unknown')
    folder_path = os.path.join("Results", model_name)
    os.makedirs(folder_path, exist_ok=True)
    return folder_path


def _save_sim_data_to_csv(sim_data, folder_path=None, filename='sim_data.csv'):
    """
    Save simulation data to a CSV file.

    Args:
        sim_data: The simulation data object, typically from MPC interface.
        folder_path: Directory where the CSV file will be saved. Defaults to current working directory.
        filename: Name of the CSV file. Defaults to 'sim_data.csv'.

    Returns:
        None
    """
    sim_data_export = sim_data.get_data()
    if isinstance(sim_data, mpc.data.scalar_data.ScalarData):
        sim_time = {'Time': ['0']}
    else:
        sim_time = {'Time': sim_data.get_time_points() or []}
    data = {**sim_time, **sim_data_export}

    folder_path = folder_path or os.getcwd()
    os.makedirs(folder_path, exist_ok=True)

    file_path = os.path.join(folder_path, filename)
    try:
        with open(file_path, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(data.keys())

            if isinstance(sim_data, mpc.data.scalar_data.ScalarData):
                row = [data[key] for key in data.keys()]
                writer.writerow(row)
            else:
                max_length = max(len(value) for value in data.values())
                for i in range(max_length):
                    row = [data[key][i] if i < len(data[key]) else '' for key in data.keys()]
                    writer.writerow(row)
        print(f"[Info] Simulation data saved to {file_path}")
    except Exception as e:
        print(f"[Error] Could not save simulation data to {file_path}: {e}")


def _save_figure(figure, folder_path=None, filename='figure.png', dpi=300):
    """
    Save a matplotlib figure to the specified folder.

    Args:
        figure: Matplotlib figure object to save.
        folder_path: Directory path to save the figure.
        filename: Name of the output file.
        dpi: Dots-per-inch for raster formats. Defaults to 300.

    Returns:
        None
    """
    folder_path = folder_path or os.getcwd()
    os.makedirs(folder_path, exist_ok=True)
    file_path = os.path.join(folder_path, filename)
    try:
        _, file_extension = os.path.splitext(filename)
        file_extension = file_extension.lower()
        vector_formats = ['.pdf', '.svg', '.eps']
        if file_extension in vector_formats:
            figure.savefig(file_path)
        else:
            figure.savefig(file_path, dpi=dpi)
        print(f"[Info] Figure saved to {file_path}")
    except Exception as e:
        print(f"[Error] Could not save figure to {file_path}: {e}")


def _save_io_csv(time_series, io_data_array, plant, folder_path, filename='io_data.csv'):
    """
    Write the full CV/MV trajectory to CSV, overwriting on each call.

    Column headers are prefixed with ``CV_`` or ``MV_`` (using the variable's
    display name) so that the ``live_plot.py`` watcher can identify them
    without any extra metadata.

    Parameters
    ----------
    time_series : list of float
        Simulation time at each row.
    io_data_array : list of list
        Each entry is ``[cv_values..., mv_values...]`` for one time step.
    plant : Plant
        Provides ``CV_index``, ``MV_index``, ``CV_display_names``, and
        ``MV_display_names``.
    folder_path : str
        Directory where the file is written.
    filename : str, optional
        Output filename. Defaults to ``'io_data.csv'``.
    """
    cv_headers = [f"CV_{plant.CV_display_names[j]}" for j in range(len(plant.CV_index))]
    mv_headers = [f"MV_{plant.MV_display_names[j]}" for j in range(len(plant.MV_index))]
    headers = ["Time"] + cv_headers + mv_headers

    os.makedirs(folder_path, exist_ok=True)
    file_path = os.path.join(folder_path, filename)
    try:
        with open(file_path, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(headers)
            for t, row in zip(time_series, io_data_array):
                writer.writerow([t] + ['' if v is None else v for v in row])
    except Exception as e:
        print(f"[Error] Could not save IO data to {file_path}: {e}")


def _plot_final_results(time_series, io_data_array, plant, show=False):
    """
    Generate final CV and MV trajectory plots.

    Args:
        time_series: List of time points.
        io_data_array: 2D list of simulation results (CVs + MVs).
        plant: The plant model.
        show: Whether to immediately display the plots.

    Returns:
        figures: List of matplotlib figure objects.
        names: Corresponding sanitized figure name identifiers.
    """
    io_np = np.array(io_data_array)
    time_label = f"${plant.time_display_name[0]}$"
    figures = []
    names = []

    for j, var_name in enumerate(plant.CV_index):
        fig = plt.figure()
        display_name = plant.CV_display_names[j]
        plt.plot(time_series, io_np[:, j], label="Trajectory")
        sp = plant.steady_state_values.get(var_name, None)
        if sp is not None:
            plt.axhline(y=sp, color='r', linestyle='--', label="Setpoint")
        plt.ylabel(f"${display_name}$")
        plt.xlabel(time_label)
        plt.title(f"$\\mathrm{{CV}}_{{{j+1}}}:\\;{display_name}\\;\\mathrm{{vs.}}\\;\\mathrm{{Time}}$")
        plt.grid(True)
        plt.legend()
        figures.append(fig)
        names.append(f"CV_{j+1}_{display_name}")

    for j, var_name in enumerate(plant.MV_index):
        fig = plt.figure()
        display_name = plant.MV_display_names[j]
        mv_index = len(plant.CV_index) + j
        plt.plot(time_series, io_np[:, mv_index], label="Trajectory")
        sp = plant.steady_state_values.get(var_name, None)
        if sp is not None:
            plt.axhline(y=sp, color='r', linestyle='--', label="Setpoint")
        plt.ylabel(f"${display_name}$")
        plt.xlabel(time_label)
        plt.title(f"$\\mathrm{{MV}}_{{{j+1}}}:\\;{display_name}\\;\\mathrm{{vs.}}\\;\\mathrm{{Time}}$")
        plt.grid(True)
        plt.legend()
        figures.append(fig)
        names.append(f"MV_{j+1}_{display_name}")

    return figures, names


def _handle_mpc_results(sim_data, time_series, io_data_array, plant, cpu_time, options,
                        folder_path=None):
    """
    Post-processing after the MPC loop: saves data, plots, and manages output directories.

    Parameters
    ----------
    sim_data : TimeSeriesData
        Full plant trajectory.
    time_series : list of float
        Sampling-interval time stamps.
    io_data_array : list of list
        CV + MV values at each sampling instant.
    plant : Plant
        Provides variable metadata for labelling.
    cpu_time : list of float
        Solver CPU time per iteration.
    options : Options
        Simulation configuration.
    folder_path : str or None
        Directory to write results into.  If None, falls back to a legacy
        path based on options (no timestamp).
    """
    final_figures, figure_names = _plot_final_results(time_series, io_data_array, plant)

    if folder_path is None:
        # Fallback: reconstruct a legacy path (no timestamp).
        folder_path = _get_results_folder(options)

    try:
        os.makedirs(folder_path, exist_ok=True)
    except Exception as e:
        print(f"[Error] Could not create results directory: {folder_path}\n{e}")
        return

    if options.save_data:
        _save_sim_data_to_csv(sim_data, folder_path=folder_path)
        _save_io_csv(time_series, io_data_array, plant, folder_path)
        _write_plot_script(folder_path, options, plant, sim_data)

    if options.save_figure:
        for idx, fig in enumerate(final_figures):
            raw_name = figure_names[idx] if idx < len(figure_names) else f"figure_{idx+1}"
            name = re.sub(r'[^a-zA-Z0-9_\-]', '_', raw_name)
            _save_figure(fig, folder_path=folder_path, filename=f"{name}.png")

    plt.close('all')

    average_cpu_time = sum(cpu_time) / len(cpu_time) if cpu_time else 0
    print(f"Average CPU time per solve: {average_cpu_time:.4f} seconds")
    print(f"Results saved to: {folder_path}")


import csv

def _save_epsilon(iteration, LHS, options):
    """
    Append the minimum epsilon (LHS) value for one MPC iteration to a CSV file.

    Creates ``LHS_values.csv`` in the results folder on the first call and
    appends subsequent rows. The LHS value is the minimum value of epsilon in
    [0, 1) that certifies closed-loop stability for that horizon.

    Parameters
    ----------
    iteration : int
        Zero-based MPC iteration index.
    LHS : float
        The computed LHS / minimum epsilon value for this iteration.
    options : Options
        Simulation configuration used to locate the results folder.
    """
    folder_path = _get_results_folder(options)
    csv_filename = os.path.join(folder_path, "LHS_values.csv")

    # Create the file and header if it doesn't exist
    if not os.path.exists(csv_filename):
        with open(csv_filename, mode="w", newline="") as file:
            writer = csv.writer(file)
            writer.writerow(["iteration", "LHS"])

    # Append new LHS value
    with open(csv_filename, mode="a", newline="") as file:
        writer = csv.writer(file)
        writer.writerow([iteration, LHS])

    
def _plot_lyap(lyap, options):
    """
    Plot the Lyapunov function value trajectory over the MPC horizon.

    Displays a blocking figure of the Lyapunov value at each horizon step,
    useful for verifying closed-loop stability during post-run analysis.

    Parameters
    ----------
    lyap : list of float
        Lyapunov function values recorded at each MPC iteration.
    options : Options
        Simulation configuration (currently unused but included for consistency).
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import threading

    lyap = np.asarray(lyap)[1:]
    x = np.arange(1, len(lyap) + 1)

    print("Plotting Lyapunov Function")

    plt.figure()
    plt.plot(x, lyap)
    plt.xlabel("Number of Horizons")
    plt.ylabel("Lyapunov Function Value")
    plt.tight_layout()
    plt.show(block=True)


def _save_lyap_csv(lyap, options, filename="lyapunov.csv"):
    """
    Save the Lyapunov function trajectory to a CSV file.

    Parameters
    ----------
    lyap : list of float
        Lyapunov function values recorded at each MPC iteration.
    options : Options
        Simulation configuration used to locate the results folder.
    filename : str, optional
        Output filename. Defaults to ``'lyapunov.csv'``.
    """
    import numpy as np
    import pandas as pd
    import os

    folder_path = _get_results_folder(options)
    full_path = os.path.join(folder_path, filename)

    lyap = np.asarray(lyap)[1:]

    df = pd.DataFrame({
        "Horizon": np.arange(1, len(lyap) + 1),
        "Lyapunov Value": lyap
    })

    df.to_csv(full_path, index=False)
