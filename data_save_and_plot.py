import os
import csv
import pyomo.contrib.mpc as mpc
import matplotlib.pyplot as plt
import numpy as np
import re


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

    average_cpu_time = sum(cpu_time) / len(cpu_time) if cpu_time else 0
    print(f"Average CPU time per solve: {average_cpu_time:.4f} seconds")

    if options.plot_end:
        plt.show()