import pandas as pd
import matplotlib.pyplot as plt
import os
import glob
from matplotlib.ticker import MaxNLocator

# ========================================
# Custom Settings
# ========================================
data_dir = 'Plotting Tools'
custom_run_labels = [r"$N=2$", r"$N=20$", r"$N=2/\infty$"]

figures = [
    {
        "label": "Fig1",
        "header_pattern": "Ca[*]",
        "xlabel": "Time (min)",
        "ylabel": r"$C_a\;(mol/L)$",
        "output_filename": "Summary_Fig1.png",
        "include_setpoint": True,
        "setpoint_value": 0.5,
        "setpoint_label": "Steady State Optimum"
    },
    {
        "label": "Fig2",
        "header_pattern": "Cb[*]",
        "xlabel": "Time (min)",
        "ylabel": r"$C_b\;(mol/L)$",
        "output_filename": "Summary_Fig2.png",
        "include_setpoint": True,
        "setpoint_value": 0.5,
        "setpoint_label": "Steady State Optimum"
    },
    {
        "label": "Fig3",
        "header_pattern": "Fa0[*]",
        "xlabel": "Time (min)",
        "ylabel": r"$F_{A0}\;(mol/min)$",
        "output_filename": "Summary_Fig3.png",
        "include_setpoint": True,
        "setpoint_value": 12,
        "setpoint_label": "Steady State Optimum"
    }
]

def match_columns(columns, pattern):
    return [col for col in columns if col == pattern]

# Load and plot each figure
csv_files = sorted(glob.glob(os.path.join(data_dir, '*.csv')))

for fig in figures:
    fig_obj, ax = plt.subplots(figsize=(7, 5))
    run_idx = 0
    time_data = None  # store time for setpoint plotting

    for file in csv_files:
        df = pd.read_csv(file)
        df.columns = [col.strip() for col in df.columns]

        time_col = next((col for col in df.columns if 'time' in col.lower()), None)
        if time_col is None:
            continue

        matching_cols = match_columns(df.columns, fig["header_pattern"])
        if not matching_cols:
            continue

        if time_data is None:
            time_data = df[time_col]  # use first valid file to get time axis

        for col in matching_cols:
            label = custom_run_labels[run_idx] if run_idx < len(custom_run_labels) else f"Run {run_idx + 1}"
            linestyle = '--' if label == r"$N=2/\infty$" else '-'
            ax.plot(df[time_col], df[col], label=label, linestyle=linestyle, linewidth=2.5)
            run_idx += 1

    # Plot custom setpoint line
    if fig.get("include_setpoint") and time_data is not None:
        setpoint_label = fig.get("setpoint_label", "Setpoint")
        ax.plot(time_data, [fig["setpoint_value"]] * len(time_data),
                '--', color='black', linewidth=2.5, label=setpoint_label)

    if run_idx == 0 and not fig.get("include_setpoint"):
        print(f"⚠️  No data plotted for {fig['label']}!")
        continue

    # Set axis labels and ticks
    ax.set_xlabel(fig["xlabel"], fontsize=16)
    ax.set_ylabel(fig["ylabel"], fontsize=16)
    ax.tick_params(axis='both', which='major', labelsize=14)

    # Limit number of ticks
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))

    ax.legend(fontsize=14, frameon=False)
    fig_obj.tight_layout()

    # Save figure at high resolution
    fig_obj.savefig(os.path.join(data_dir, fig["output_filename"]), dpi=1200, bbox_inches='tight')

plt.show()
