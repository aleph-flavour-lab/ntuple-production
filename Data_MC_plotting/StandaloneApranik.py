import uproot
import awkward as ak
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
import run_list

# ======================================================================
#  1. LATEX & PLOT CONFIGURATION
# ======================================================================
def configure_plots():
    tex_bin = "/cvmfs/sft.cern.ch/lcg/external/texlive/2020/bin/x86_64-linux"
    if os.path.isdir(tex_bin):
        os.environ["PATH"] = tex_bin + ":" + os.environ.get("PATH", "")
    has_latex = os.system("which latex > /dev/null 2>&1") == 0
    rcParams.update({
        "text.usetex": has_latex, 
        "font.family": "serif", 
        "font.size": 10,
        "axes.labelsize": 11, 
        "figure.figsize": (5, 5), 
        "figure.autolayout": True
    })

configure_plots()

# Physics Constants
LUMI, XSEC_ZQQ, N_GEN_TOT = run_list.luminosity_pb(1994), 30385.0, 1080986  # pb^-1 of the runs kept by stage1
WEIGHT = (XSEC_ZQQ * LUMI) / N_GEN_TOT

DATA_DIR = "/afs/cern.ch/work/h/hfatehi/Birgit/stage1/data/"
MC_DIR   = "/afs/cern.ch/work/h/hfatehi/Birgit/stage1/"
OUTPUT_DIR = "/eos/user/h/hfatehi/aleph_final_plots2"

FLAVORS = [
    ("Zdd", r"$Z \to d\bar{d}$", "#D6E4F0"),
    ("Zuu", r"$Z \to u\bar{u}$", "#85C1E9"),
    ("Zss", r"$Z \to s\bar{s}$", "#2E86C1"),
    ("Zcc", r"$Z \to c\bar{c}$", "#1A5276"),
    ("Zbb", r"$Z \to b\bar{b}$", "#0F407B"),
]

# ======================================================================
#  2. FULL HARD-CODED VARIABLE MAP
# ======================================================================
# Format: "var_name": (xmin, xmax, bins, label)
VAR_MAP = {
    # --- Jet Variables ---
    "jet_mass":   (0, 60, 60, r"Jet Mass [GeV]"),
    "jet_p":      (10, 80, 70, r"Jet Momentum $p$ [GeV]"),
    "jet_e":      (10, 80, 70, r"Jet Energy $E$ [GeV]"),
    "jet_phi":    (-3.14, 3.14, 50, r"Jet $\phi$ [rad]"),
    "jet_theta":  (0.6, 2.6, 50, r"Jet $\theta$ [rad]"),
    "jet_pT":     (0, 80, 80, r"Jet $p_{T}$ [GeV]"),
    "jet_nnhad":  (0, 17, 17, r"Number of Neutral Hadrons"),
    "jet_ngamma": (0, 25, 25, r"Number of Photons"),
    "jet_nchad":  (0, 30, 30, r"Number of Charged Hadrons"),
    "jet_nel":    (0, 6, 6,   r"Number of Electrons"),
    "jet_nmu":    (0, 5, 5,   r"Number of Muons"),
    "jet_nconst": (0, 55, 55, r"Jet Particle Multiplicity"),

    # --- PFCand PID & Booleans ---
    "pfcand_isMu": (0, 2, 2, r"isMuon Flag"),
    "pfcand_isEl": (0, 2, 2, r"isElectron Flag"),
    "pfcand_isChargedHad": (0, 2, 2, r"isChargedHad Flag"),
    "pfcand_isGamma": (0, 2, 2, r"isGamma Flag"),
    "pfcand_isNeutralHad": (0, 2, 2, r"isNeutralHad Flag"),

    # --- PFCand Kinematics ---
    "pfcand_e": (0, 50, 50, r"Jet Constituents Energy $E$ [GeV]"),
    "pfcand_p": (0, 50, 50, r"Jet Constituents Momentum $p$ [GeV]"),
    "pfcand_theta": (0, 3.14, 30, r"Jet Constituents $\theta$ [rad]"),
    "pfcand_phi": (-3.14, 3.14, 30, r"Jet Constituents $\phi$ [rad]"),
    "pfcand_charge": (-1.5, 1.5, 3, r"Jet Constituents Charge"),
    #"pfcand_type": (0, 10, 10, r"Jet Constituents Type"),
    "pfcand_erel": (0, 1, 30, r"Relative Energy $E_{cand}/E_{jet}$"),
    "pfcand_erel_log": (-2.6, 0, 30, r"$\log(E_{rel})$"),
    "pfcand_thetarel": (0, 0.5, 30, r"Relative $\theta$ to Jet"),
    "pfcand_phirel": (0, 0.5, 30, r"Relative $\phi$ to Jet"),

    # --- Impact Parameters ---
    "pfcand_dxy": (-0.5, 0.5, 30, r"Impact Parameter $d_{xy}$ "),
    "pfcand_dz": (-0.5, 0.5, 30, r"Impact Parameter $d_{z}$ "),
    "pfcand_phi0": (-3.14, 3.14, 30, r"Track $\phi_0$"),
    "pfcand_C": (-0.005, 0.025, 30, r"Track Curvature $C$"),
    "pfcand_ct": (-5, 5, 30, r"Track $cot(\theta)$"),

    # --- Covariances ---
    "pfcand_dptdpt": (0, 0.002, 30, r"$cov(p_T, p_T)$"),
    "pfcand_dxydxy": (0, 0.001, 30, r"$cov(d_{xy}, d_{xy})$"),
    "pfcand_dzdz": (0, 0.005, 30, r"$cov(d_z, d_z)$"),
    "pfcand_dphidphi": (0, 0.0001, 30, r"$cov(\phi, \phi)$"),
    "pfcand_detadeta": (0, 0.0001, 30, r"$cov(\eta, \eta)$"),

    # --- B-Tagging ---
    "pfcand_btagSip2dVal": (-0.1, 0.5, 30, r"2D IP Value "),
    "pfcand_btagSip2dSig": (-30, 30, 30, r"2D IP Significance"),
    "pfcand_btagSip3dVal": (-0.1, 1.0, 30, r"3D IP Value "),
    "pfcand_btagSip3dSig": (-30, 30, 30, r"3D IP Significance"),
    "pfcand_btagJetDistVal": (0, 0.1, 30, r"Distance to Jet Axis "),
    "pfcand_btagJetDistSig": (-20, 20, 30, r"Distance to Jet Axis Significance"),

    # --- dE/dx & PID (Pads) ---
    "pfcand_dEdx_pads_value": (0, 20, 30, r"dE/dx (Pads) [[GeV]]"),
    "pfcand_PID_pval_pads_ele": (-1, 1, 30, r"Electron $p$-val (Pads)"),
    "pfcand_PID_pval_pads_mu": (-1, 1, 30, r"Muon $p$-val (Pads)"),
    "pfcand_PID_pval_pads_pi": (-1, 1, 30, r"Pion $p$-val (Pads)"),
    "pfcand_PID_pval_pads_kaon": (-1, 1, 30, r"Kaon $p$-val (Pads)"),
    "pfcand_PID_pval_pads_proton": (-1, 1, 30, r"Proton $p$-val (Pads)"),

    # --- dE/dx & PID (Wires) ---
    "pfcand_dEdx_wires_value": (0, 15, 30, r"dE/dx (Wires) [GeV]"),
    "pfcand_PID_pval_wires_ele": (-1, 1, 30, r"Electron $p$-val (Wires)"),
    "pfcand_PID_pval_wires_mu": (-1, 1, 30, r"Muon $p$-val (Wires)"),
    "pfcand_PID_pval_wires_pi": (-1, 1, 30, r"Pion $p$-val (Wires)"),
    "pfcand_PID_pval_wires_kaon": (-1, 1, 30, r"Kaon $p$-val (Wires)"),
    "pfcand_PID_pval_wires_proton": (-1, 1, 30, r"Proton $p$-val (Wires)"),
}

# ======================================================================
#  3. PLOTTING ENGINE
# ======================================================================
def get_flat_array(file_path, var_name):
    if not os.path.exists(file_path): return None
    with uproot.open(file_path + ":events") as tree:
        data = tree.arrays([var_name, "jet_p", "jet_theta"], library="ak")
        mask = (ak.all(data["jet_p"] > 10.0, axis=1) & 
                ak.all(abs(np.cos(data["jet_theta"])) < 0.65, axis=1))
        return ak.flatten(data[var_name][mask], axis=None).to_numpy()

def plot_variable(var_name):
    if var_name not in VAR_MAP:
        print(f"Skipping {var_name} (not in map)")
        return
        
    xmin, xmax, bins, label = VAR_MAP[var_name]
    data_vals = get_flat_array(os.path.join(DATA_DIR, "1994.root"), var_name)
    if data_vals is None or len(data_vals) == 0: return

    mc_hists, mc_colors, mc_labels, mc_weights, all_mc_flat = [], [], [], [], []
    for fname, name, color in FLAVORS:
        vals = get_flat_array(os.path.join(MC_DIR, f"{fname}.root"), var_name)
        if vals is not None:
            mc_hists.append(vals)
            mc_colors.append(color)
            mc_labels.append(name)
            mc_weights.append(np.ones_like(vals) * WEIGHT)
            all_mc_flat.append(vals)

    d_counts, edges = np.histogram(data_vals, bins=bins, range=(xmin, xmax))
    centers = (edges[:-1] + edges[1:]) / 2
    mc_total_vals = np.concatenate(all_mc_flat)
    mc_counts, _ = np.histogram(mc_total_vals, bins=bins, range=(xmin, xmax), weights=np.ones_like(mc_total_vals) * WEIGHT)

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, height_ratios=[3, 1])
    plt.subplots_adjust(hspace=0.07)

    # Labels outside at height 1.01
    ax1.text(0.0, 1.01, r"\textbf{ALEPH} Archived Data", transform=ax1.transAxes, ha='left', va='bottom')
    ax1.text(1.0, 1.01, r"$\sqrt{s} = 91.2$ GeV, \ $\mathcal{L} = %.2f \ \mathrm{pb}^{-1}$" % LUMI, transform=ax1.transAxes, ha='right', va='bottom')

    # Main Plot
    ax1.hist(mc_hists, bins=bins, range=(xmin, xmax), weights=mc_weights, 
             stacked=True, color=mc_colors, label=mc_labels, alpha=0.8, histtype='stepfilled')
    ax1.errorbar(centers, d_counts, yerr=np.sqrt(d_counts), fmt='ko', markersize=2, label="Data - 1994")

    ax1.set_yscale('log')
    ax1.set_ylim(1e2, 1e8)
    ax1.set_ylabel("Events")
    ax1.legend(ncol=2, fontsize=7, loc='upper right', frameon=False)

    # Ratio Plot
    ratio = np.divide(d_counts, mc_counts, out=np.zeros_like(d_counts, dtype=float), where=mc_counts!=0)
    ax2.errorbar(centers, ratio, yerr=np.divide(np.sqrt(d_counts), mc_counts, out=np.zeros_like(ratio), where=mc_counts!=0), fmt='ko', markersize=2)
    ax2.axhline(1.0, color='red', lw=0.8, ls='--')
    ax2.set_ylabel("Data/MC")
    ax2.set_xlabel(label)
    ax2.set_ylim(0.5, 1.5)
    ax2.grid(axis='y', alpha=0.2)

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    plt.savefig(f"{OUTPUT_DIR}/{var_name}.png", bbox_inches='tight', dpi=600)
    plt.close()

if __name__ == "__main__":
    for v in VAR_MAP.keys():
        print(f"Plotting {v}...")
        try:
            plot_variable(v)
        except Exception as e:
            print(f"Error on {v}: {e}")
