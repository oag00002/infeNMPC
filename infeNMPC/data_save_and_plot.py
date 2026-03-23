import os
import csv
import pyomo.contrib.mpc as mpc
import matplotlib.pyplot as plt
import numpy as np
import re


def _get_results_folder(options):
    """
    Construct and create the hierarchical results directory for the current run.

    The path encodes the key simulation parameters so that results from
    different configurations are automatically segregated:
    ``Results/infinite_horizon_{true|false}/finite_horizon_{N}/
    gamma_{g}/beta_{b}/disturbance_{true|false}/``.

    Parameters
    ----------
    options : Options
        Simulation configuration used to determine the folder path.

    Returns
    -------
    str
        Absolute (or CWD-relative) path to the results directory.
    """
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


def _handle_mpc_results(sim_data, time_series, io_data_array, plant, cpu_time, options):
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

    plt.close('all')

    average_cpu_time = sum(cpu_time) / len(cpu_time) if cpu_time else 0
    print(f"Average CPU time per solve: {average_cpu_time:.4f} seconds")


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
