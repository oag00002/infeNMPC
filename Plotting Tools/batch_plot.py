import pandas as pd
import matplotlib.pyplot as plt
import os
import glob
from matplotlib.ticker import MaxNLocator

# ========================================
# Custom Settings
# ========================================
data_dir = 'Plotting Tools/Distillation eNMPC Results 3'
custom_run_labels = [r"$N=3,tracking$", r"$N=20,tracking$", r"$N-1=2,tracking$",
                     r"$N=3,economic$", r"$N=20,economic$", r"$N-1=2,economic$"]

figures = [
    {
        "label": "Fig1",
        "header_pattern": "x[1,*]",
        "xlabel": "Time (min)",
        "ylabel": r"$(1-x_1)$",
        "output_filename": "Summary_Fig1.png",
        "include_setpoint": True,
        "setpoint_label": "Setpoint"
    },
    {
        "label": "Fig2",
        "header_pattern": "T[29,*]",
        "xlabel": "Time (min)",
        "ylabel": r"$T_{29}\;(K)$",
        "output_filename": "Summary_Fig2.png",
        "include_setpoint": True,
        "setpoint_label": "Setpoint"
    },
    {
        "label": "Fig3",
        "header_pattern": "Qr[*]",
        "xlabel": "Time (min)",
        "ylabel": r"$Reboiler\;Heat\;Duty\;(MJ)$",
        "output_filename": "Summary_Fig3.png",
        "include_setpoint": False
    },
    {
        "label": "Fig4",
        "header_pattern": "Rec[*]",
        "xlabel": "Time (min)",
        "ylabel": r"$Reflux\;Ratio$",
        "output_filename": "Summary_Fig4.png",
        "include_setpoint": False
    }
]

def match_columns(columns, pattern):
    return [col for col in columns if col == pattern]

# Load and plot each figure
csv_files = sorted(glob.glob(os.path.join(data_dir, '*.csv')))

for fig in figures:
    fig_obj, ax = plt.subplots(figsize=(7, 5))  # larger figure for better resolution
    run_idx = 0

    for file in csv_files:
        df = pd.read_csv(file)
        df.columns = [col.strip() for col in df.columns]

        time_col = next((col for col in df.columns if 'time' in col.lower()), None)
        if time_col is None:
            continue

        matching_cols = match_columns(df.columns, fig["header_pattern"])
        if not matching_cols:
            continue

        if fig.get("include_setpoint") and "Setpoint" in df.columns and run_idx == 0:
            ax.plot(df[time_col], df["Setpoint"], '--', color='black', linewidth=2.5,
                    label=fig["setpoint_label"])

        for col in matching_cols:
            label = custom_run_labels[run_idx] if run_idx < len(custom_run_labels) else f"Run {run_idx + 1}"
            linestyle = '--' if label == r"$N-1=2$" else '-'
            ax.plot(df[time_col], df[col], label=label, linestyle=linestyle, linewidth=2.5)
            run_idx += 1

    if run_idx == 0:
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
