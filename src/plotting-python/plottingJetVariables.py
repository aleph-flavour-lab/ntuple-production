import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import os
import run_list

# ============================================
# --- ROOT / HEP-like Matplotlib Style ---
# ============================================
plt.rcParams.update({
    "text.usetex": False,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "axes.labelsize": 18,
    "axes.titlesize": 20,
    "axes.linewidth": 1.2,
    "axes.labelweight": "normal",
    "xtick.labelsize": 15,
    "ytick.labelsize": 15,
    "legend.fontsize": 14,
    "figure.dpi": 300,

    # ROOT-like markers & lines
    "lines.linewidth": 2.0,
    "lines.markersize": 6,

    # ROOT-like grid
    "grid.color": "0.85",
    "grid.linestyle": "--",
    "grid.linewidth": 0.7,

    "errorbar.capsize": 2,
    "patch.edgecolor": "black",
    "patch.linewidth": 0.8,
})

# ROOT-like MC palette
# colors = [
#     "#1a85ff",  # blue
#     "#28a745",  # green
#     "#ff7f0e",  # orange
#     "#984ea3",  # violet
#     "#c0392b",  # red/brown
# ]
#colors = ['#1f77b4', '#3fa9a9', '#5dc0c0', '#80d1d1', '#a3e0e0']
colors = ['#a3e0e0', '#80d1d1', '#5dc0c0', '#3fa9a9', '#1f77b4']

# ============================================
# --- Files and MC info ---
# ============================================
files = {
    r"$Z \to d\bar{d}$": {"path": "/eos/user/h/hfatehi/D0fliped/Stage1/Zdd.root", "type": "mc"},
    r"$Z \to u\bar{u}$": {"path": "/eos/user/h/hfatehi/D0fliped/Stage1/Zuu.root", "type": "mc"},
    r"$Z \to s\bar{s}$": {"path": "/eos/user/h/hfatehi/D0fliped/Stage1/Zss.root", "type": "mc"},
    r"$Z \to c\bar{c}$": {"path": "/eos/user/h/hfatehi/D0fliped/Stage1/Zcc.root", "type": "mc"},
    r"$Z \to b\bar{b}$": {"path": "/eos/user/h/hfatehi/D0fliped/Stage1/Zbb.root", "type": "mc"},
    "LEP1-1994 (Class 16–17)": {"path": "/eos/user/h/hfatehi/D0fliped/Stage1/Znn.root", "type": "data"},
}

L_data = run_list.luminosity_pb(1994)  # pb^-1, runs kept by stage1
mc_xsec = {
    r"$Z \to b\bar{b}$": 2.81e5,
    r"$Z \to c\bar{c}$": 2.81e5 * 0.8,
    r"$Z \to s\bar{s}$": 2.81e5,
    r"$Z \to d\bar{d}$": 2.81e5,
    r"$Z \to u\bar{u}$": 2.81e5 * 0.8,
}

# ============================================
# --- Plotting Function ---
# ============================================
def plot_variable_stacked_custom(var, bins=50, logy=False, save=False, xmin=None, xmax=None, density=False):

    print(f"📊 Plotting (custom stacked): {var}")

    mc_data = {}
    data_arr = None

    # -------- Load data ----------
    for label, info in files.items():
        with uproot.open(info["path"]) as f:
            arr = f["events"][var].array(library="ak")

        if isinstance(arr, ak.Array) and arr.ndim > 1:
            arr = ak.flatten(arr)

        arr = np.array(arr, dtype=float)
        arr = arr[np.isfinite(arr)]

        if info["type"] == "data":
            data_arr = arr
        else:
            mc_data[label] = arr

    if data_arr is None:
        raise RuntimeError("❌ No data sample found.")

    # -------- Histogram Range --------
    if xmin is None:
        xmin = np.min(data_arr)
    if xmax is None:
        xmax = np.max(data_arr)

    edges = np.linspace(xmin, xmax, bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])

    # -------- MC stacking order --------
    stack_order = [
        r"$Z \to d\bar{d}$",
        r"$Z \to u\bar{u}$",
        r"$Z \to s\bar{s}$",
        r"$Z \to c\bar{c}$",
        r"$Z \to b\bar{b}$"
    ]

    mc_counts = []
    mc_labels = []
    mc_colors = []
    mc_weights = []

    for i, label in enumerate(stack_order):
        arr = mc_data[label]
        sigma = mc_xsec[label]
        N_gen = len(arr)
        w = (L_data * sigma) / N_gen

        mc_weights.append(np.ones_like(arr) * w)
        mc_counts.append(arr)
        mc_labels.append(label)
        mc_colors.append(colors[i % len(colors)])

    # -------- Global MC scaling --------
    total_mc_raw = sum(np.sum(w) for w in mc_weights)
    scale_factor = len(data_arr) / total_mc_raw
    mc_weights = [w * scale_factor for w in mc_weights]

    # =========================================
    # --- Figure + Axes ---
    # =========================================
    fig = plt.figure(figsize=(14, 14))
    gs = GridSpec(2, 1, height_ratios=[3, 1], hspace=0.05)
    ax = fig.add_subplot(gs[0])
    ax_ratio = fig.add_subplot(gs[1], sharex=ax)

    # -------- Stacked MC --------
    ax.hist(
        mc_counts, bins=edges, weights=mc_weights,
        stacked=True, color=mc_colors, label=mc_labels,
        edgecolor="black"
    )

    # -------- Data points --------
    y_data, _ = np.histogram(data_arr, bins=edges)
    y_err = np.sqrt(y_data)

    ax.errorbar(
        centers, y_data, yerr=y_err,
        fmt="o", color="black", markersize=5,
        label="LEP1-1994 (Class 16–17)"
    )

    # -------- Stats box --------
    stats = []
    stats.append(f"Data: N={len(data_arr):,}  μ={np.mean(data_arr):.3f}  σ={np.std(data_arr):.3f}")

    for i, arr in enumerate(mc_counts):
        stats.append(f"{mc_labels[i]}: μ={np.mean(arr):.3f}  σ={np.std(arr):.3f}")

    ax.text(
        0.5, 1.02, "   |   ".join(stats),
        transform=ax.transAxes, fontsize=10,
        ha="center", va="bottom",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="white", alpha=0.85)
    )

    ax.set_ylabel("Entries", fontweight="bold")
    ax.grid(alpha=0.4)
    ax.legend(frameon=False)

    if logy:
        ax.set_yscale("log")

    # =========================================
    # --- Ratio plot ---
    # =========================================
    mc_sum, _ = np.histogram(
        np.concatenate(mc_counts),
        bins=edges,
        weights=np.concatenate(mc_weights)
    )

    mask = mc_sum > 0
    ratio = np.zeros_like(mc_sum, dtype=float)
    ratio_err = np.zeros_like(mc_sum, dtype=float)
    ratio[mask] = y_data[mask] / mc_sum[mask]
    ratio_err[mask] = y_err[mask] / mc_sum[mask]

    ax_ratio.errorbar(
        centers, ratio, yerr=ratio_err,
        fmt="o", color="black", markersize=5
    )
    ax_ratio.axhline(1.0, color="red", linestyle="--", linewidth=1.2)

    ax_ratio.set_ylabel("Data/MC", fontweight="bold")
    ax_ratio.set_xlabel(var.replace("_", " ").title(), fontweight="bold")
    ax_ratio.grid(alpha=0.4)
    ax_ratio.set_ylim(0, 2)

    # =========================================
    # --- FIXED SAVE PATH ---
    # =========================================
    if save:
        eos_dir = "/eos/user/h/hfatehi/D0fliped/Stage1/jets/"
        os.makedirs(eos_dir, exist_ok=True)

        out_path = f"{eos_dir}{var}_PRD_stacked_custom.png"
        plt.savefig(out_path, bbox_inches="tight")
        print(f"✅ Saved to: {out_path}")

    plt.show()


# ============================================
# Run Plots
# ============================================
plot_variable_stacked_custom("event_invariant_mass", bins=50, xmin=50, xmax=120, density=True, save=True)
plot_variable_stacked_custom("jet_e", bins=90, xmin=0, xmax=90, density=True, save=True)
plot_variable_stacked_custom("jet_p", bins=90, xmin=0, xmax=90, density=True, save=True)
plot_variable_stacked_custom("jet_pT", bins=90, xmin=0, xmax=90, density=True, save=True)
plot_variable_stacked_custom("jet_ngamma", bins=25, xmin=0, xmax=25, density=True, save=True)
plot_variable_stacked_custom("jet_nnhad", bins=15, xmin=0, xmax=15, density=True, save=True)
plot_variable_stacked_custom("jet_nchad", bins=25, xmin=0, xmax=25, density=True, save=True)
plot_variable_stacked_custom("jet_nel", bins=5, xmin=0, xmax=5, density=True, save=True)
plot_variable_stacked_custom("jet_nmu", bins=5, xmin=0, xmax=5, density=True, save=True)
plot_variable_stacked_custom("jet_eta", bins=50, xmin=-5, xmax=5, density=True, save=True)
plot_variable_stacked_custom("jet_phi", bins=50, xmin=-10, xmax=10, density=True, save=True)
plot_variable_stacked_custom("jet_theta", bins=50, xmin=-4, xmax=4, density=True, save=True)
