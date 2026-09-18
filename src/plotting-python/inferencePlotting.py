import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import run_list

# Set Physical Review D style with LaTeX
plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Times New Roman', 'Computer Modern Roman'],
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 11,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.dpi': 300,
    'savefig.dpi': 600,
    'text.usetex': True,  # Enable LaTeX
    'text.latex.preamble': r'\usepackage{amsmath}',
    'axes.linewidth': 0.8,
    'xtick.major.width': 0.8,
    'ytick.major.width': 0.8,
    'xtick.minor.width': 0.6,
    'ytick.minor.width': 0.6,
    'lines.linewidth': 1.0,
})

mc_file = "/eos/user/h/hfatehi/alephdec9/MC-Dec08.root"
data_file = "/eos/user/h/hfatehi/alephdec9/Data-Dec08.root"

# Flavor configuration: branch_name -> label
flavors = {
    "recojet_isD": r"$D$-jets",
    "recojet_isU": r"$U$-jets",
    "recojet_isS": r"$S$-jets",
    "recojet_isC": r"$C$-jets",
    "recojet_isB": r"$B$-jets",
}
# PRD-appropriate color scheme (colorblind-friendly)
colors = ['#8dd3c7', '#bebada', '#fb8072', '#80b1d3', '#fdb462']

# Cross-section normalization factor
xs_normalization = run_list.luminosity_pb(1994) * 30385 / 1080986  # pb^-1 of the runs kept by stage1
print(f"Cross-section normalization factor: {xs_normalization:.6f}")

# ------------------------
# Load MC B-scores with theta filter
# ------------------------
with uproot.open(mc_file) as f:
    bscore_mc = f["Events"]["score_recojet_isB"].array(library="ak")
    jet_theta_mc = f["Events"]["jet_theta"].array(library="ak")
    
    # Flatten if needed
    if bscore_mc.ndim > 1:
        bscore_mc = ak.flatten(bscore_mc)
    if jet_theta_mc.ndim > 1:
        jet_theta_mc = ak.flatten(jet_theta_mc)
    
    bscore_mc = ak.to_numpy(bscore_mc)
    jet_theta_mc = ak.to_numpy(jet_theta_mc)
    
    # Apply theta cut: 0.75 < |theta| < 2.4
    theta_mask_mc = (np.abs(jet_theta_mc) > 0.75) & (np.abs(jet_theta_mc) < 2.4)
    bscore_mc = bscore_mc[theta_mask_mc]

mc_counts = []
mc_labels = []

with uproot.open(mc_file) as f:
    for flavor_branch, flavor_label in flavors.items():
        # Load flavor truth
        is_flavor = f["Events"][flavor_branch].array(library="ak")
        jet_theta = f["Events"]["jet_theta"].array(library="ak")
        
        if is_flavor.ndim > 1:
            is_flavor = ak.flatten(is_flavor)
        if jet_theta.ndim > 1:
            jet_theta = ak.flatten(jet_theta)
        
        is_flavor = ak.to_numpy(is_flavor).astype(bool)
        jet_theta = ak.to_numpy(jet_theta)
        
        # Apply theta cut
        theta_mask = (np.abs(jet_theta) > 0.75) & (np.abs(jet_theta) < 2.4)
        combined_mask = is_flavor & theta_mask
        
        # Select jets of this flavor with theta cut
        scores_flavor = bscore_mc[combined_mask[theta_mask_mc]]
        scores_flavor = scores_flavor[np.isfinite(scores_flavor)]
        
        # Ensure even number of jets (for dijets)
        n_jets = len(scores_flavor)
        n_pairs = n_jets // 2
        scores_flavor = scores_flavor[:n_pairs*2]
        
        # Compute dijet scores: score_jet1 * score_jet2
        dijet_scores = scores_flavor[0::2] * scores_flavor[1::2]
        mc_counts.append(dijet_scores)
        mc_labels.append(flavor_label)
        print(f"{flavor_label}: {n_pairs} dijets, mean dijet score = {np.mean(dijet_scores):.3f}")

# ------------------------
# Load Data B-scores with theta filter
# ------------------------
with uproot.open(data_file) as f:
    bscore_data = f["Events"]["score_recojet_isB"].array(library="ak")
    jet_theta_data = f["Events"]["jet_theta"].array(library="ak")
    
    if bscore_data.ndim > 1:
        bscore_data = ak.flatten(bscore_data)
    if jet_theta_data.ndim > 1:
        jet_theta_data = ak.flatten(jet_theta_data)
    
    bscore_data = ak.to_numpy(bscore_data)
    jet_theta_data = ak.to_numpy(jet_theta_data)
    
    # Apply theta cut: 0.75 < |theta| < 2.4
    theta_mask_data = (np.abs(jet_theta_data) > 0.75) & (np.abs(jet_theta_data) < 2.4)
    bscore_data = bscore_data[theta_mask_data]
    bscore_data = bscore_data[np.isfinite(bscore_data)]

# Compute dijet scores for data
n_jets_data = len(bscore_data)
n_pairs_data = n_jets_data // 2
bscore_data = bscore_data[:n_pairs_data*2]
dijet_data = bscore_data[0::2] * bscore_data[1::2]
print("Data: N dijets =", len(dijet_data), "mean dijet score =", np.mean(dijet_data))

# ------------------------
# Compute logit transformation
# ------------------------
def safe_logit(x, epsilon=1e-7):
    """Compute logit with clipping to avoid infinities"""
    x_clipped = np.clip(x, epsilon, 1 - epsilon)
    return np.log(x_clipped / (1 - x_clipped))

dijet_data_logit = safe_logit(dijet_data)
mc_counts_logit = [safe_logit(arr) for arr in mc_counts]

# ------------------------
# Histogram bins for both representations
# ------------------------
# Original score bins
bins_score = 20
edges_score = np.linspace(0, 1, bins_score + 1)
centers_score = 0.5 * (edges_score[:-1] + edges_score[1:])

# Logit bins
bins_logit = 20
logit_min, logit_max = -8, 8
edges_logit = np.linspace(logit_min, logit_max, bins_logit + 1)
centers_logit = 0.5 * (edges_logit[:-1] + edges_logit[1:])

# ------------------------
# Apply cross-section normalization to MC
# ------------------------
mc_weights = [np.ones_like(arr) * xs_normalization for arr in mc_counts]
mc_weights_logit = [np.ones_like(arr) * xs_normalization for arr in mc_counts_logit]

# Print normalization info
mc_total_before = sum(len(arr) for arr in mc_counts)
mc_total_after = mc_total_before * xs_normalization
print(f"\nMC total events before normalization: {mc_total_before}")
print(f"MC total events after XS normalization: {mc_total_after:.2f}")
print(f"Data total events: {len(dijet_data)}")
print(f"MC/Data ratio: {mc_total_after/len(dijet_data):.3f}")

# ========================
# PLOT 1: Original Score
# ========================
def compute_ratio_with_errors(data_hist, mc_hist):
    """Compute ratio with Poisson errors"""
    ratio = np.zeros_like(data_hist, dtype=float)
    ratio_err = np.zeros_like(data_hist, dtype=float)
    
    for i in range(len(data_hist)):
        if mc_hist[i] > 0:
            ratio[i] = data_hist[i] / mc_hist[i]
            if data_hist[i] > 0:
                ratio_err[i] = ratio[i] * np.sqrt(data_hist[i]) / data_hist[i]
        else:
            ratio[i] = 0
            ratio_err[i] = 0
    
    return ratio, ratio_err

# Compute histograms for score
y_data_score, _ = np.histogram(dijet_data, bins=edges_score)
y_mc_total_score = np.zeros(bins_score)
for arr, wgt in zip(mc_counts, mc_weights):
    h, _ = np.histogram(arr, bins=edges_score, weights=wgt)
    y_mc_total_score += h

ratio_score, ratio_err_score = compute_ratio_with_errors(y_data_score, y_mc_total_score)

# Create figure for score
fig1 = plt.figure(figsize=(3.4, 4.2))  # PRD single column width
gs1 = GridSpec(2, 1, height_ratios=[3, 1], hspace=0.08)

# Main plot - Score
ax1 = fig1.add_subplot(gs1[0])
ax1.hist(
    mc_counts,
    bins=edges_score,
    weights=mc_weights,
    stacked=True,
    color=colors,
    label=mc_labels,
    linewidth=0,
    edgecolor='none'
)

y_err_score = np.sqrt(y_data_score)
ax1.errorbar(
    centers_score, y_data_score, yerr=y_err_score,
    fmt='o', color='black', markersize=3, linewidth=0.8, 
    label='Data', zorder=10, capsize=0
)

ax1.set_ylabel(r"Events", fontsize=11)
ax1.set_yscale('log')
ax1.set_ylim(1, 2e10)
ax1.set_xlim(0, 1)
ax1.legend(loc='upper right', fontsize=8, frameon=False)
ax1.tick_params(axis='x', labelbottom=False, direction='in', top=True, right=True)
ax1.tick_params(axis='y', direction='in', right=True)
ax1.minorticks_on()

# Add ALEPH label with LaTeX
ax1.text(
    0.05, 0.95,
    r'\textbf{ALEPH} Preliminary',
    transform=ax1.transAxes,
    fontsize=10,
    verticalalignment='top'
)

ax1.text(
    0.05, 0.88,
    (
        r'$0.75 < |\theta_{\mathrm{jet}}| < 2.4$' "\n"
        r'$p_T^{\mathrm{Jet}} > 10\ \mathrm{GeV}$'
    ),
    transform=ax1.transAxes,
    fontsize=9,
    verticalalignment='top'
)


# Ratio plot - Score
ax2 = fig1.add_subplot(gs1[1], sharex=ax1)
ax2.errorbar(
    centers_score, ratio_score, yerr=ratio_err_score,
    fmt='o', color='black', markersize=3, linewidth=0.8, capsize=0
)
ax2.axhline(y=1, color='red', linestyle='--', linewidth=0.8)
ax2.fill_between([0, 1], 0.8, 1.2, color='yellow', alpha=0.3)

ax2.set_xlabel(r"Dijet $b$-tag score ($\mathrm{score}_{j_1} \times \mathrm{score}_{j_2}$)", fontsize=11)
ax2.set_ylabel(r"Data/MC", fontsize=10)
ax2.set_ylim(0, 2)
ax2.set_xlim(0, 1)
ax2.grid(True, alpha=0.2, linewidth=0.5, linestyle=':')
ax2.tick_params(direction='in', top=True, right=True)
ax2.minorticks_on()

plt.savefig("/eos/user/h/hfatehi/prd_dijet_score.pdf", bbox_inches='tight', pad_inches=0.02)
plt.savefig("/eos/user/h/hfatehi/prd_dijet_score.png", bbox_inches='tight', pad_inches=0.02)

plt.show()

print("\nPlot saved:")
print("  - prd_dijet_score.pdf/png (original score with theta cut)")
