import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import os
import gc
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import run_list

# ================================================================
#  PAPER-STYLE PLOTTING (ROOT-like but cleaner)
# ================================================================
plt.rcParams.update({
    "text.usetex": False,
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
    "axes.labelweight": "bold",
    "axes.titleweight": "bold",
    "axes.labelsize": 17,
    "axes.titlesize": 19,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "legend.fontsize": 13,
    "figure.dpi": 300,
    "axes.linewidth": 1.4,
    "lines.linewidth": 1.8,
    "errorbar.capsize": 0,
    "grid.color": "0.85",
    "grid.linewidth": 1.0,
})

# Soft, paper-friendly stacked colors
colors = [
    "#c7e9f1", "#9ed9de", "#72c7c7",
    "#4bb0ae", "#2a8b9b"
]

# ================================================================
#  FILES & MC INFO
# ================================================================
eos = "/eos/user/h/hfatehi/D0fliped/Stage1/"

files = {
    "LEP1-1994 (Class 16–17)": {"path": eos + "Znn.root", "type": "data"},
    r"$Z \to d\bar{d}$": {"path": eos + "Zdd.root", "type": "mc"},
    r"$Z \to u\bar{u}$": {"path": eos + "Zuu.root", "type": "mc"},
    r"$Z \to s\bar{s}$": {"path": eos + "Zss.root", "type": "mc"},
    r"$Z \to c\bar{c}$": {"path": eos + "Zcc.root", "type": "mc"},
    r"$Z \to b\bar{b}$": {"path": eos + "Zbb.root", "type": "mc"},
}

stack_order = [
    r"$Z \to d\bar{d}$",
    r"$Z \to u\bar{u}$",
    r"$Z \to s\bar{s}$",
    r"$Z \to c\bar{c}$",
    r"$Z \to b\bar{b}$"
]

# Lumi & x-sections
L_data = run_list.luminosity_pb(1994)  # pb^-1, runs kept by stage1
mc_xsec = {
    r"$Z \to b\bar{b}$": 2.81e5,
    r"$Z \to c\bar{c}$": 2.81e5 * 0.8,
    r"$Z \to s\bar{s}$": 2.81e5,
    r"$Z \to d\bar{d}$": 2.81e5,
    r"$Z \to u\bar{u}$": 2.81e5 * 0.8,
}

# ================================================================
#  STACKED PLOT FUNCTION (OPTIMIZED)
# ================================================================
def plot_variable_stacked(var, bins=50, logy=False, save=False, xmin=None, xmax=None):

    print(f"\n📊 Plotting variable: {var}")

    mc_data = {}
    data_arr = None

    # ---------------------------
    # Load ROOT files variable-wise
    # ---------------------------
    for label, info in files.items():

        with uproot.open(info["path"]) as f:
            arr = f["events"][var].array(library="ak")

        arr = ak.to_numpy(ak.flatten(arr, axis=None))
        arr = arr[np.isfinite(arr)]  # safety

        if info["type"] == "data":
            data_arr = arr
        else:
            mc_data[label] = arr

    if data_arr is None:
        raise RuntimeError("❌ ERROR: No data array loaded!")

    # ---------------------------
    # Histogram edges
    # ---------------------------
    if xmin is None: xmin = np.min(data_arr)
    if xmax is None: xmax = np.max(data_arr)

    edges = np.linspace(xmin, xmax, bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])

    # ---------------------------
    # MC weights
    # ---------------------------
    mc_counts = []
    mc_weights = []
    mc_colors = []

    for i, label in enumerate(stack_order):
        arr = mc_data[label]
        sigma = mc_xsec[label]
        N = len(arr)
        w = (L_data * sigma) / N

        mc_counts.append(arr)
        mc_weights.append(np.full_like(arr, w, dtype=float))
        mc_colors.append(colors[i])

    # Global MC → Data scaling
    total_mc_events = sum(np.sum(w) for w in mc_weights)
    scale_factor = len(data_arr) / total_mc_events
    mc_weights = [w * scale_factor for w in mc_weights]

    # ---------------------------
    # Figure layout
    # ---------------------------
    fig = plt.figure(figsize=(13, 13))
    gs = GridSpec(2, 1, height_ratios=[3, 1], hspace=0.03)

    ax = fig.add_subplot(gs[0])
    ax_ratio = fig.add_subplot(gs[1], sharex=ax)

    # ---------------------------
    # MC stacked
    # ---------------------------
    ax.hist(
        mc_counts,
        bins=edges,
        weights=mc_weights,
        stacked=True,
        color=mc_colors,
        edgecolor="black",
        linewidth=0.7,
        label=stack_order,
    )

    # ---------------------------
    # DATA
    # ---------------------------
    y_data, _ = np.histogram(data_arr, bins=edges)
    y_err = np.sqrt(y_data)

    ax.errorbar(
        centers, y_data, yerr=y_err,
        fmt="o", color="black", markersize=3.5,
        label="LEP1-1994 (Class 16–17)"
    )

    # ---------------------------
    # RATIO
    # ---------------------------
    mc_total, _ = np.histogram(
        np.concatenate(mc_counts),
        bins=edges,
        weights=np.concatenate(mc_weights)
    )

    ratio = np.divide(
        y_data, mc_total,
        out=np.zeros_like(y_data, dtype=float),
        where=mc_total > 0
    )
    ratio_err = np.divide(
        y_err, mc_total,
        out=np.zeros_like(y_err, dtype=float),
        where=mc_total > 0
    )

    ax_ratio.errorbar(
        centers, ratio, yerr=ratio_err,
        fmt="o", color="black", markersize=3
    )
    ax_ratio.axhline(1.0, color="red", linestyle="--", linewidth=1.2)

    ax_ratio.set_ylim(0.4, 1.6)

    # ---------------------------
    # Labels & aesthetics
    # ---------------------------
    ax.set_ylabel("Entries", fontweight="bold")
    ax_ratio.set_ylabel("Data/MC", fontweight="bold")
    ax_ratio.set_xlabel(var.replace("_", " ").capitalize(), fontweight="bold")

    ax.grid(alpha=0.4)
    ax_ratio.grid(alpha=0.4)

    if logy:
        ax.set_yscale("log")

    plt.setp(ax.get_xticklabels(), visible=False)
    ax.legend(frameon=False, loc="upper right")

    # ---------------------------
    # Save
    # ---------------------------
    if save:
        outdir = eos + "pfcand/"
        os.makedirs(outdir, exist_ok=True)
        outname = f"{outdir}/{var}_stacked.png"
        plt.savefig(outname, bbox_inches="tight")
        print(f"✅ Saved to: {outname}")

    plt.show()
    plt.close(fig)

    print(f"📌 Data events: {len(data_arr):,}")
    print(f"📌 Scaled MC sum: {np.sum(mc_total):,.2f}")

    # ---------------------------
    # MEMORY CLEANUP
    # ---------------------------
    del mc_data, data_arr, mc_counts, mc_weights
    del y_data, y_err, mc_total, ratio, ratio_err
    gc.collect()


# ================================================================
#  RUN PLOTS
# ================================================================
plot_variable_stacked("pfcand_dxy", bins=100, logy=True, save=True, xmin=-1, xmax=1)
plot_variable_stacked("pfcand_dz",  bins=100, logy=True, save=True, xmin=-3, xmax=3)
