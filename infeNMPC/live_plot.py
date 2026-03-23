"""
Live plotter for infeNMPC safe_run mode.

Watches an ``io_data.csv`` written by the MPC loop (when ``safe_run=True``)
and redraws CV/MV trajectories whenever the file is updated.  Columns whose
headers start with ``CV_`` are plotted in the left column; ``MV_`` columns go
in the right column.

Usage
-----
Run in a separate terminal while the MPC loop is running::

    python -m infeNMPC.live_plot <path/to/io_data.csv>

The plotter polls every second and updates automatically.  Close the window
or press Ctrl-C to stop.
"""
import subprocess
import sys
import os
import time
import csv

import matplotlib.pyplot as plt


def start_live_plotter(options):
    """
    Launch the live plotter as a background subprocess and return the handle.

    Computes the ``io_data.csv`` path from *options* (the same location that
    ``safe_run`` writes to), then spawns ``python -m infeNMPC.live_plot`` in a
    separate process.  The caller is responsible for calling
    ``proc.terminate()`` / ``proc.wait()`` when the MPC run finishes.

    Parameters
    ----------
    options : Options
        Simulation configuration.  Must have ``safe_run=True`` for the CSV to
        be written during the run.

    Returns
    -------
    subprocess.Popen
        Handle to the plotter process.

    Example
    -------
    ::

        options = Options.for_model_module(model, safe_run=True, ...)
        plotter = start_live_plotter(options)
        try:
            mpc_loop(options)
        finally:
            plotter.terminate()
            plotter.wait()
    """
    from .data_save_and_plot import _get_results_folder
    csv_path = os.path.join(_get_results_folder(options), 'io_data.csv')
    proc = subprocess.Popen(
        [sys.executable, '-m', 'infeNMPC.live_plot', csv_path],
    )
    return proc


def _read_csv(path):
    with open(path, newline='') as f:
        rows = list(csv.DictReader(f))
    if not rows:
        return [], {}, {}
    headers = list(rows[0].keys())
    cv_cols = [h for h in headers if h.startswith('CV_')]
    mv_cols = [h for h in headers if h.startswith('MV_')]
    times = [float(r['Time']) for r in rows]
    cv_data = {
        h: [float(r[h]) if r[h] != '' else None for r in rows]
        for h in cv_cols
    }
    mv_data = {
        h: [float(r[h]) if r[h] != '' else None for r in rows]
        for h in mv_cols
    }
    return times, cv_data, mv_data


def show_final_results(options):
    """
    Read the saved ``io_data.csv`` and display a blocking final plot.

    Call this after ``mpc_loop`` finishes (and after terminating the live
    plotter subprocess) to show the complete CV/MV trajectories in the same
    two-column layout as the live monitor.  The window stays open until the
    user closes it.

    Parameters
    ----------
    options : Options
        Simulation configuration — used to locate the results folder.
    """
    from .data_save_and_plot import _get_results_folder
    csv_path = os.path.join(_get_results_folder(options), 'io_data.csv')

    if not os.path.exists(csv_path):
        print(f"[Warning] No io_data.csv found at {csv_path}; skipping final plot.")
        return

    times, cv_data, mv_data = _read_csv(csv_path)
    num_cv = len(cv_data)
    num_mv = len(mv_data)
    num_rows = max(num_cv, num_mv, 1)

    fig, axes = plt.subplots(
        num_rows, 2,
        figsize=(14, 4 * num_rows),
        squeeze=False,
        sharex=True,
    )
    fig.suptitle("MPC Final Results")

    for j, (col, vals) in enumerate(cv_data.items()):
        ax = axes[j, 0]
        ax.plot(times, vals)
        ax.set_title(col[3:])
        ax.set_ylabel(col[3:])
        ax.set_xlabel("Time")
        ax.grid(True)
    for j in range(num_cv, num_rows):
        axes[j, 0].set_visible(False)

    for j, (col, vals) in enumerate(mv_data.items()):
        ax = axes[j, 1]
        ax.plot(times, vals)
        ax.set_title(col[3:])
        ax.set_ylabel(col[3:])
        ax.set_xlabel("Time")
        ax.grid(True)
    for j in range(num_mv, num_rows):
        axes[j, 1].set_visible(False)

    plt.tight_layout()
    plt.show(block=True)


def main():
    if len(sys.argv) < 2:
        print("Usage: python -m infeNMPC.live_plot <path/to/io_data.csv>")
        sys.exit(1)

    csv_path = sys.argv[1]

    print(f"Watching {csv_path} ...")
    while not os.path.exists(csv_path):
        time.sleep(0.5)

    times, cv_data, mv_data = _read_csv(csv_path)
    num_cv = len(cv_data)
    num_mv = len(mv_data)
    num_rows = max(num_cv, num_mv, 1)

    plt.ion()
    fig, axes = plt.subplots(
        num_rows, 2,
        figsize=(14, 4 * num_rows),
        squeeze=False,
        sharex=True,
    )
    fig.suptitle("infeNMPC Live Monitor")

    last_mtime = 0.0

    def _redraw(times, cv_data, mv_data):
        for j, (col, vals) in enumerate(cv_data.items()):
            ax = axes[j, 0]
            ax.cla()
            ax.plot(times, vals)
            ax.set_title(col[3:])   # strip "CV_"
            ax.set_ylabel(col[3:])
            ax.set_xlabel("Time")
            ax.grid(True)
        for j in range(len(cv_data), num_rows):
            axes[j, 0].set_visible(False)

        for j, (col, vals) in enumerate(mv_data.items()):
            ax = axes[j, 1]
            ax.cla()
            ax.plot(times, vals)
            ax.set_title(col[3:])   # strip "MV_"
            ax.set_ylabel(col[3:])
            ax.set_xlabel("Time")
            ax.grid(True)
        for j in range(len(mv_data), num_rows):
            axes[j, 1].set_visible(False)

        plt.tight_layout()
        fig.canvas.draw()
        fig.canvas.flush_events()

    _redraw(times, cv_data, mv_data)

    try:
        while plt.get_fignums():
            try:
                mtime = os.path.getmtime(csv_path)
            except OSError:
                plt.pause(1.0)
                continue
            if mtime > last_mtime:
                last_mtime = mtime
                try:
                    times, cv_data, mv_data = _read_csv(csv_path)
                    _redraw(times, cv_data, mv_data)
                except Exception as e:
                    print(f"[Warning] Could not read {csv_path}: {e}")
            plt.pause(1.0)
    except KeyboardInterrupt:
        pass


if __name__ == '__main__':
    main()
