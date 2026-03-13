import os
import csv
import io
import sys
import pyomo.contrib.mpc as mpc
import matplotlib.pyplot as plt
import numpy as np
import re


def _get_results_folder(options):
    if options.infinite_horizon:
        folder_path = os.path.join(
            "Results",
            "infinite_horizon_true",
            f"finite_horizon_{options.nfe_finite}",
            f"gamma_{options.gamma}",
            f"beta_{options.beta}",
            f"disturbance_{str(options.disturb_flag).lower()}"
        )
    else:
        folder_path = os.path.join(
            "Results",
            "infinite_horizon_false",
            f"finite_horizon_{options.nfe_finite}",
            f"disturbance_{str(options.disturb_flag).lower()}"
        )

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


def _setup_live_plot(plant):
    """
    Set up a live plot for visualizing MV and CV variables.

    Args:
        plant: The plant model with MV and CV definitions.

    Returns:
        fig: The matplotlib figure.
        axes: The axes array for subplots.
    """
    plt.ion()
    num_MV = len(plant.MV_display_names)
    num_CV = len(plant.CV_display_names)
    num_rows = max(num_MV, num_CV)

    screen_width = 1920
    screen_height = 1080
    dpi = 100
    fig_width = screen_width / dpi
    fig_height = screen_height / dpi

    fig, axes = plt.subplots(num_rows, 2, figsize=(fig_width, fig_height), sharex=True)
    fig.suptitle("Live MPC Simulation")
    return fig, axes


def _update_live_plot(fig, axes, time_series, io_data_array, plant):
    """
    Update the live plot with current simulation data.

    Args:
        fig: Matplotlib figure object.
        axes: Array of subplot axes.
        time_series: List of simulation time points.
        io_data_array: List of CV and MV data at each time.
        plant: The plant model with display names and steady state values.

    Returns:
        None
    """
    io_np = np.array(io_data_array)
    time_label = f"${plant.time_display_name[0]}$"

    # Ensure axes is 2D for uniform indexing
    axes = np.atleast_2d(axes)

    for j, var_name in enumerate(plant.CV_index):
        ax = axes[j, 0] if axes.shape[1] > 1 else axes[j]
        ax.cla()
        display_name = plant.CV_display_names[j]
        ax.plot(time_series, io_np[:, j], label="Trajectory")
        sp = plant.steady_state_values.get(var_name, None)
        if sp is not None:
            ax.axhline(y=sp, color='r', linestyle='--', label="Setpoint")
        ax.set_ylabel(f"${display_name}$")
        ax.set_xlabel(time_label)
        ax.set_title(f"$\\mathrm{{CV}}_{{{j+1}}}:\\;{display_name}\\;\\mathrm{{vs.}}\\;\\mathrm{{Time}}$")
        ax.grid(True)
        ax.legend()

    for j, var_name in enumerate(plant.MV_index):
        ax = axes[j, 1] if axes.shape[1] > 1 else axes[j + len(plant.CV_index)]
        ax.cla()
        display_name = plant.MV_display_names[j]
        mv_index = len(plant.CV_index) + j
        ax.plot(time_series, io_np[:, mv_index], label="Trajectory")
        sp = plant.steady_state_values.get(var_name, None)
        if sp is not None:
            ax.axhline(y=sp, color='r', linestyle='--', label="Setpoint")
        ax.set_ylabel(f"${display_name}$")
        ax.set_xlabel(time_label)
        ax.set_title(f"$\\mathrm{{MV}}_{{{j+1}}}:\\;{display_name}\\;\\mathrm{{vs.}}\\;\\mathrm{{Time}}$")
        ax.grid(True)
        ax.legend()

    plt.tight_layout()
    plt.pause(0.01)


def _finalize_live_plot(fig):
    """
    Finalize and close the interactive live plot.

    Args:
        fig: Matplotlib figure object.

    Returns:
        None
    """
    plt.ioff()
    plt.close(fig)


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


def _handle_mpc_results(sim_data, time_series, io_data_array, plant, solver_time, options):
    """
    Post-processing after the MPC loop: saves data, plots, and manages output directories.

    Args:
        sim_data: Simulation data object.
        time_series: List of time steps.
        io_data_array: Combined CV and MV trajectories.
        plant: The plant model.
        options: Configuration object with output flags and parameters.

    Returns:
        None
    """
    final_figures, figure_names = _plot_final_results(time_series, io_data_array, plant)

    if options.infinite_horizon:
        folder_path = os.path.join(
            "Results",
            "infinite_horizon_true",
            f"finite_horizon_{options.nfe_finite}",
            f"gamma_{options.gamma}",
            f"beta_{options.beta}",
            f"disturbance_{str(options.disturb_flag).lower()}"
        )
    else:
        folder_path = os.path.join(
            "Results",
            "infinite_horizon_false",
            f"finite_horizon_{options.nfe_finite}",
            f"disturbance_{str(options.disturb_flag).lower()}"
        )
    try:
        os.makedirs(folder_path, exist_ok=True)
    except Exception as e:
        print(f"[Error] Could not create results directory: {folder_path}\n{e}")
        return

    if options.save_data:
        _save_sim_data_to_csv(sim_data, folder_path=folder_path)

    if options.save_figure:
        for idx, fig in enumerate(final_figures):
            raw_name = figure_names[idx] if idx < len(figure_names) else f"figure_{idx+1}"
            name = re.sub(r'[^a-zA-Z0-9_\-]', '_', raw_name)
            _save_figure(fig, folder_path=folder_path, filename=f"{name}.png")

    average_solver_time = sum(solver_time) / len(solver_time) if solver_time else 0
    print(f"Average CPU time per solve: {average_solver_time:.4f} seconds")

    if options.plot_end:
        plt.show()


import csv

def _save_epsilon(iteration, LHS, options):
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


# ---------------------------------------------------------------------------
# IPOPT iteration log capture and saving
# ---------------------------------------------------------------------------

_IPOPT_COLUMNS = ["iter", "objective", "inf_pr", "inf_du", "lg_mu", "norm_d", "lg_rg", "alpha_du", "alpha_pr", "ls"]
# Matches a single IPOPT iteration line: integer then exactly 9 more whitespace-separated tokens.
_IPOPT_ITER_RE = re.compile(
    r'^\s*(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$'
)


def _parse_ipopt_log(log_text):
    """Parse captured IPOPT output and return a list of dicts, one per iteration."""
    rows = []
    for line in log_text.splitlines():
        m = _IPOPT_ITER_RE.match(line)
        if m:
            # IPOPT appends a step-type letter to alpha_pr (e.g. '1.09e-02f'); strip it.
            alpha_pr = re.sub(r'[a-zA-Z]+$', '', m.group(9))
            rows.append({
                "iter":      int(m.group(1)),
                "objective": m.group(2),
                "inf_pr":    m.group(3),
                "inf_du":    m.group(4),
                "lg_mu":     m.group(5),
                "norm_d":    m.group(6),
                "lg_rg":     m.group(7),
                "alpha_du":  m.group(8),
                "alpha_pr":  alpha_pr,
                "ls":        m.group(10).strip(),
            })
    return rows


def _solve_and_log_ipopt(solver, model, mpc_iter, label, tee_flag, log_dir="ipopt_logs"):
    """
    Solve *model* with *solver*, capture the IPOPT iteration log, and write it
    to ``<log_dir>/iter_<mpc_iter:04d>_<label>.csv``.

    Parameters
    ----------
    solver    : Pyomo SolverFactory instance
    model     : Pyomo model to solve
    mpc_iter  : integer MPC loop index (used for file naming)
    label     : string identifier, e.g. 'controller' or 'plant'
    tee_flag  : bool — if True, solver output is also printed to the console
    log_dir   : directory for CSV files (created if absent)

    Returns
    -------
    Pyomo solver result object
    """
    buf = io.StringIO()

    class _Tee:
        """Mirror writes to both the real stdout and a capture buffer."""
        def __init__(self, stream, buf):
            self._stream = stream
            self._buf = buf
        def write(self, data):
            self._stream.write(data)
            self._buf.write(data)
        def flush(self):
            self._stream.flush()
            self._buf.flush()
        def __getattr__(self, name):
            return getattr(self._stream, name)

    old_stdout = sys.stdout
    sys.stdout = _Tee(old_stdout, buf) if tee_flag else buf
    try:
        result = solver.solve(model, tee=True)
    finally:
        sys.stdout = old_stdout

    rows = _parse_ipopt_log(buf.getvalue())

    os.makedirs(log_dir, exist_ok=True)
    path = os.path.join(log_dir, f"iter_{mpc_iter:04d}_{label}.csv")
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=_IPOPT_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)

    return result
