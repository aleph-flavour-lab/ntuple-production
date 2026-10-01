#!/usr/bin/env python3
"""
ALEPH 1994 Z -> qqbar: overplot three MC productions per flavour.

Productions: Paper (legacy) | HVFL | Pythia 8.3 (Monash)

For each flavour in order (Zbb, Zcc, Zss, Zuu, Zdd):
  * read the three productions' ROOT files
  * fill the same set of histograms
  * normalise each histogram to unit area (AUC = 1)
  * overplot the three normalized distributions + ratio panels

No weights, no N_gen, no data, no luminosity.
"""

import argparse
import concurrent.futures
import glob
import json
import os
import time

import awkward as ak
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams
import numpy as np
import uproot


# ======================================================================
# 1. LATEX & PLOT CONFIGURATION
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
        "figure.figsize": (6, 5),
        "figure.autolayout": False,
    })
    return has_latex


HAS_LATEX = configure_plots()


# ======================================================================
# 2. CONFIG
# ======================================================================

NUM_WORKERS = 64
STEP_SIZE = "100 MB"
DPI = 300
SCALES = ("linear", "log")

# (key, legend label, base directory, colour).  First entry = ratio reference.
PRODUCTIONS = [
    ("Paper", r"Paper (legacy)",
     "/eos/experiment/fcc/ee/analyses/case-studies/aleph/"
     "processedMC/1994/zqq/stage1/comparison/legacyMC",
     "#1A5276"),
    ("HVFL", r"HVFL",
     "/eos/experiment/fcc/ee/analyses/case-studies/aleph/"
     "processedMC/1994/zqq/stage1/comparison/hvfl-40M",
     "#C0392B"),
    ("Pythia8.3", r"Pythia 8.3 (Monash)",
     "/eos/experiment/fcc/ee/analyses/case-studies/aleph/"
     "processedMC/1994/zqq/stage1/comparison/v-pythia8-1.0",
     "#1E8449"),
]

# Order in which flavours are processed (and grouped on disk).
FLAVORS = ["Zbb", "Zcc", "Zss", "Zuu", "Zdd"]

FLAVOR_LABELS = {
    "Zbb": r"$Z \to b\bar{b}$",
    "Zcc": r"$Z \to c\bar{c}$",
    "Zss": r"$Z \to s\bar{s}$",
    "Zuu": r"$Z \to u\bar{u}$",
    "Zdd": r"$Z \to d\bar{d}$",
}

OUTPUT_BASE = "./mc_compare"
RANGE_CACHE = "./var_map_suggested_mc_compare.py"


# ======================================================================
# 3. PRESELECTION  (unchanged semantics: event- and jet-level, per branch)
# ======================================================================

OBJECT_PREFIXES = {
    "constituent": "pfcand_",
    "sv": "sv_",
    "v0": "v0_",
}

PRESELECTION = {
    "event": {
        "event_njet": {"min": 1, "label": r"$N_{\rm jet} \geq 2$"},
    },
    "jet": {
        "jet_pT": {"min": 10.0, "require": "all",
                   "label": r"$p_{T,\rm jet} > 10$ GeV"},
        "jet_eta": {"max": 0.9, "transform": "abs", "require": "all",
                    "label": r"$|\eta_{\rm jet}| < 0.9$"},
    },
    "constituent": {},
    "sv": {},
    "v0": {},
}


# ======================================================================
# 4. VARIABLE MAP  -  (xmin, xmax, nbins, axis label)
# ======================================================================

VAR_MAP = {
    # ------------------------------------------------- event variables -----
    "event_class": (0, 32, 32, r"Event class bit index"),
    "event_number": (0, 200000, 100, r"Event number"),
    "run_number": (14000, 32000, 90, r"Run number"),
    "event_invariant_mass": (0, 120, 60, r"$m_{jj}$ [GeV]"),
    "event_njet": (0, 6, 6, r"Number of jets"),
    "VertexX": (-1.0, 1.0, 50, r"Stored vertex $x$"),
    "VertexY": (-1.0, 1.0, 50, r"Stored vertex $y$"),
    "VertexZ": (-10.0, 10.0, 50, r"Stored vertex $z$"),

    # ------------------------------------------- primary vertex refit ------
    "n_primary_tracks": (0, 40, 40, r"$N_{\rm tracks}^{\rm primary}$"),
    "n_secondary_tracks": (0, 30, 30, r"$N_{\rm tracks}^{\rm secondary}$"),
    "Beamspot_x": (-100, 100, 50, r"Beamspot $x$ [$10\,\mu$m]"),
    "Beamspot_y": (-100, 100, 50, r"Beamspot $y$ [$10\,\mu$m]"),
    "Beamspot_z": (-100, 100, 50, r"Beamspot $z$"),
    "Beamspot_x_cm": (-0.1, 0.1, 50, r"Beamspot $x$ [cm]"),
    "Beamspot_y_cm": (-0.1, 0.1, 50, r"Beamspot $y$ [cm]"),
    "Beamspot_z_cm": (-0.1, 0.1, 50, r"Beamspot $z$ [cm]"),
    "Vertex_refit_x": (-1.0, 1.0, 50, r"PV$_{\rm refit}$ $x$"),
    "Vertex_refit_y": (-1.0, 1.0, 50, r"PV$_{\rm refit}$ $y$"),
    "Vertex_refit_z": (-10.0, 10.0, 50, r"PV$_{\rm refit}$ $z$"),
    "Vertex_refit_cov_xx": (0, 1e-4, 50, r"PV cov$(x,x)$ [cm$^2$]"),
    "Vertex_refit_cov_yx": (-2e-5, 2e-5, 50, r"PV cov$(y,x)$ [cm$^2$]"),
    "Vertex_refit_cov_yy": (0, 1e-4, 50, r"PV cov$(y,y)$ [cm$^2$]"),
    "Vertex_refit_cov_zx": (-2e-5, 2e-5, 50, r"PV cov$(z,x)$ [cm$^2$]"),
    "Vertex_refit_cov_zy": (-2e-5, 2e-5, 50, r"PV cov$(z,y)$ [cm$^2$]"),
    "Vertex_refit_cov_zz": (0, 1e-3, 50, r"PV cov$(z,z)$ [cm$^2$]"),
    "Vertex_refit_chi2": (0, 50, 50, r"PV fit $\chi^2$"),
    "pv_converged": (0, 2, 2, r"PV fit converged"),
    "pv_split_converged": (0, 2, 2, r"PV split fit converged"),
    "pv_trivial": (0, 2, 2, r"PV trivial (beamspot only)"),
    "pv_good": (0, 2, 2, r"PV usable (pv\_good)" if HAS_LATEX else "PV usable (pv_good)"),
    "trk_nCand": (0, 6, 6, r"Stored candidates per track"),

    # ---------------------------------- gen vertex & resolution (MC only) --
    "gen_vertex_x": (-0.5, 0.5, 50, r"PV$_{\rm gen}$ $x$"),
    "gen_vertex_y": (-0.5, 0.5, 50, r"PV$_{\rm gen}$ $y$"),
    "gen_vertex_z": (-10.0, 10.0, 50, r"PV$_{\rm gen}$ $z$"),
    "res_vertex_x": (-0.2, 0.2, 60, r"$x_{\rm PV}^{\rm reco} - x_{\rm PV}^{\rm gen}$"),
    "res_vertex_y": (-0.2, 0.2, 60, r"$y_{\rm PV}^{\rm reco} - y_{\rm PV}^{\rm gen}$"),
    "res_vertex_z": (-0.5, 0.5, 60, r"$z_{\rm PV}^{\rm reco} - z_{\rm PV}^{\rm gen}$"),

    # ------------------------------------------- secondary vertices --------
    "n_sv_event": (0, 10, 10, r"$N_{\rm SV}$ per event"),
    "n_sv_jets": (0, 8, 8, r"$N_{\rm SV}$ per jet"),
    "sv_chi2": (0, 15, 50, r"SV $\chi^2$"),
    "sv_chi2_norm": (0, 10, 50, r"SV $\chi^2/{\rm ndof}$"),
    "sv_ndof": (0, 8, 8, r"SV ndof"),
    "sv_ntracks": (2, 10, 8, r"SV track multiplicity"),
    "sv_mass": (0, 6, 60, r"SV mass [GeV]"),
    "sv_p": (0, 40, 60, r"SV momentum $p$ [GeV]"),
    "sv_thetarel": (-0.5, 0.5, 50, r"SV $\theta$ relative to jet"),
    "sv_phirel": (-0.5, 0.5, 50, r"SV $\phi$ relative to jet"),
    "sv_dxy": (0, 20, 60, r"SV $d_{xy}$ from PV"),
    "sv_dxyz": (0, 30, 60, r"SV $d_{3D}$ from PV"),
    "sv_cosPointing": (0.6, 1.0, 50, r"SV $\cos\alpha_{\rm pointing}$"),
    "sv_prel": (0, 1, 50, r"SV $p_{\rm rel}$ to jet"),
    "sv_correctedMass": (0, 10, 50, r"SV corrected mass [GeV]"),
    "sv_dx": (-3, 3, 60, r"SV $\Delta x$ from PV"),
    "sv_dy": (-3, 3, 60, r"SV $\Delta y$ from PV"),
    "sv_dz": (-3, 3, 60, r"SV $\Delta z$ from PV"),
    "sv_cov_xx": (0, 0.01, 50, r"SV cov$(x,x)$ [cm$^2$]"),
    "sv_cov_yy": (0, 0.01, 50, r"SV cov$(y,y)$ [cm$^2$]"),
    "sv_cov_zz": (0, 0.01, 50, r"SV cov$(z,z)$ [cm$^2$]"),

    # -------------------------------------------------- V0 candidates ------
    "n_v0_event": (0, 8, 8, r"$N_{V^0}$ per event"),
    "n_v0_jets": (0, 6, 6, r"$N_{V^0}$ per jet"),
    "n_v0_ks": (0, 6, 6, r"$N_{K^0_S}$ per jet"),
    "n_v0_lambda": (0, 6, 6, r"$N_{\Lambda}$ per jet"),
    "v0_pdg": (0, 3500, 70, r"$V^0$ $|{\rm PDG}|$"),
    "v0_invM": (0, 2.0, 80, r"$V^0$ invariant mass [GeV]"),
    "v0_chi2": (0, 10, 50, r"$V^0$ $\chi^2$"),
    "v0_chi2_norm": (0, 10, 50, r"$V^0$ $\chi^2/{\rm ndof}$"),
    "v0_ndof": (0, 5, 5, r"$V^0$ ndof"),
    "v0_ntracks": (2, 6, 4, r"$V^0$ track multiplicity"),
    "v0_p": (0, 30, 60, r"$V^0$ momentum $p$ [GeV]"),
    "v0_prel": (0, 1, 50, r"$V^0$ $p_{\rm rel}$ to jet"),
    "v0_thetarel": (-0.4, 0.4, 50, r"$V^0$ $\theta$ relative to jet"),
    "v0_phirel": (-0.4, 0.4, 50, r"$V^0$ $\phi$ relative to jet"),
    "v0_dxy": (0, 60, 60, r"$V^0$ $d_{xy}$ from PV"),
    "v0_dxyz": (0, 100, 50, r"$V^0$ $d_{3D}$ from PV"),
    "v0_cosPointing": (0.98, 1.0, 50, r"$V^0$ $\cos\alpha_{\rm pointing}$"),
    "v0_correctedMass": (0, 3, 60, r"$V^0$ corrected mass [GeV]"),
    "v0_dx": (-30, 30, 60, r"$V^0$ $\Delta x$ from PV"),
    "v0_dy": (-30, 30, 60, r"$V^0$ $\Delta y$ from PV"),
    "v0_dz": (-30, 30, 60, r"$V^0$ $\Delta z$ from PV"),

    # ------------------------------------------------ track variables ------
    "n_tracks_all": (0, 60, 60, r"$N_{\rm tracks}$ (all)"),
    "n_tracks_sel": (0, 50, 50, r"$N_{\rm tracks}$ (baseline selected)"),
    "n_trackstates_sel": (0, 50, 50, r"$N_{\rm trackstates}$ (selected)"),
    "n_tracks_sel_vertexfit": (0, 50, 50, r"$N_{\rm tracks}$ (vertex fit)"),
    "chi2_tracks_all": (0, 20, 50, r"Track $\chi^2$"),
    "ndf_tracks_all": (0, 30, 30, r"Track ndof"),
    "chi2_o_ndf_tracks_all": (0, 10, 50, r"Track $\chi^2/$ndof"),

    # -------------------------------------------------- jet variables ------
    "JetClustering_d23": (0, 40, 60, r"$\sqrt{d_{23}}$ [GeV]"),
    "JetClustering_d34": (0, 25, 50, r"$\sqrt{d_{34}}$ [GeV]"),
    "jet_mass": (0, 40, 60, r"Jet Mass [GeV]"),
    "jet_p": (0, 60, 60, r"Jet Momentum $p$ [GeV]"),
    "jet_e": (0, 70, 70, r"Jet Energy $E$ [GeV]"),
    "jet_phi": (-3.15, 3.15, 50, r"Jet $\phi$ [rad]"),
    "jet_theta": (0.6, 2.6, 50, r"Jet $\theta$ [rad]"),
    "jet_pT": (0, 60, 60, r"Jet $p_T$ [GeV]"),
    "jet_eta": (-1.0, 1.0, 50, r"Jet $\eta$"),
    "jet_p_leading": (0, 60, 60, r"Leading jet $p$ [GeV]"),
    "jet_e_leading": (0, 70, 70, r"Leading jet $E$ [GeV]"),
    "jet_mass_leading": (0, 40, 60, r"Leading jet mass [GeV]"),
    "jet_phi_leading": (-3.15, 3.15, 50, r"Leading jet $\phi$ [rad]"),
    "jet_theta_leading": (0.6, 2.6, 50, r"Leading jet $\theta$ [rad]"),
    "jet_pT_leading": (0, 60, 60, r"Leading jet $p_T$ [GeV]"),
    "jet_eta_leading": (-1.0, 1.0, 50, r"Leading jet $\eta$"),
    "jet_p_subleading": (0, 60, 60, r"Subleading jet $p$ [GeV]"),
    "jet_e_subleading": (0, 70, 70, r"Subleading jet $E$ [GeV]"),
    "jet_mass_subleading": (0, 40, 60, r"Subleading jet mass [GeV]"),
    "jet_phi_subleading": (-3.15, 3.15, 50, r"Subleading jet $\phi$ [rad]"),
    "jet_theta_subleading": (0.6, 2.6, 50, r"Subleading jet $\theta$ [rad]"),
    "jet_pT_subleading": (0, 60, 60, r"Subleading jet $p_T$ [GeV]"),
    "jet_eta_subleading": (-1.0, 1.0, 50, r"Subleading jet $\eta$"),
    "jet_nnhad": (0, 20, 20, r"Number of Neutral Hadrons"),
    "jet_ngamma": (0, 30, 30, r"Number of Photons"),
    "jet_nchad": (0, 30, 30, r"Number of Charged Hadrons"),
    "jet_nel": (0, 6, 6, r"Number of Electrons"),
    "jet_nmu": (0, 6, 6, r"Number of Muons"),
    "jet_nconst": (0, 60, 60, r"Jet Particle Multiplicity"),
    "jetPID": (0, 6, 6, r"Truth jet flavour ID"),

    # ------------------------------------ PF candidate flags & kinematics --
    "pfcand_isMu": (0, 2, 2, r"isMuon Flag"),
    "pfcand_isEl": (0, 2, 2, r"isElectron Flag"),
    "pfcand_isChargedHad": (0, 2, 2, r"isChargedHad Flag"),
    "pfcand_isGamma": (0, 2, 2, r"isGamma Flag"),
    "pfcand_isNeutralHad": (0, 2, 2, r"isNeutralHad Flag"),
    "pfcand_mask": (0, 2, 2, r"Constituent mask"),
    "pfcand_e": (0, 25, 50, r"Constituent energy $E$ [GeV]"),
    "pfcand_p": (0, 25, 50, r"Constituent momentum $p$ [GeV]"),
    "pfcand_px": (-15, 15, 60, r"Constituent $p_x$ [GeV]"),
    "pfcand_py": (-15, 15, 60, r"Constituent $p_y$ [GeV]"),
    "pfcand_pz": (-40, 40, 80, r"Constituent $p_z$ [GeV]"),
    "pfcand_pt": (0, 20, 50, r"Constituent $p_T$ [GeV]"),
    "pfcand_theta": (0, 3.15, 50, r"Constituent $\theta$ [rad]"),
    "pfcand_phi": (-3.15, 3.15, 50, r"Constituent $\phi$ [rad]"),
    "pfcand_charge": (-1.5, 1.5, 3, r"Constituent charge"),
    "pfcand_erel": (0, 1, 50, r"$E_{\rm cand}/E_{\rm jet}$"),
    "pfcand_erel_log": (-6, 0, 60, r"$\log(E_{\rm rel})$"),
    "pfcand_ptrel": (0, 1, 50, r"$p_{T,\rm cand}/p_{T,\rm jet}$"),
    "pfcand_ptrel_log": (-6, 0, 60, r"$\log(p_{T,\rm rel})$"),
    "pfcand_thetarel": (0, 0.6, 50, r"Constituent $\theta$ relative to jet"),
    "pfcand_phirel": (-0.6, 0.6, 60, r"Constituent $\phi$ relative to jet"),

    # ----------------------------------------- per-constituent track info --
    "pfcand_trackChi2": (-2, 40, 42, r"Constituent track $\chi^2$"),
    "pfcand_trackNdof": (-2, 40, 42, r"Constituent track ndof"),
    "pfcand_trackChi2Norm": (-2, 10, 60, r"Constituent track $\chi^2/{\rm ndof}$"),
    "pfcand_nTrackHits_VDET": (-1, 5, 6, r"VDET hits"),
    "pfcand_nTrackHits_ITC": (-1, 9, 10, r"ITC hits"),
    "pfcand_nTrackHits_TPC": (-1, 25, 26, r"TPC hits"),

    # ------------------------------------------- impact parameters ---------
    "pfcand_dxy": (-0.5, 0.5, 50, r"Impact Parameter $d_{xy}$"),
    "pfcand_dz": (-1.0, 1.0, 50, r"Impact Parameter $d_z$"),
    "pfcand_d0": (-0.5, 0.5, 50, r"Perigee $d_0$ [cm]"),
    "pfcand_z0": (-2.0, 2.0, 50, r"Perigee $z_0$ [cm]"),
    "pfcand_phi0": (-3.15, 3.15, 50, r"Track $\phi_0$"),
    "pfcand_C": (-0.005, 0.025, 60, r"Track Curvature $C$"),
    "pfcand_ct": (-5, 5, 50, r"Track $\cot(\theta)$"),

    # ------------------------------------------------- covariances ---------
    "pfcand_dptdpt": (0, 0.002, 50, r"$\mathrm{cov}(p_T,p_T)$"),
    "pfcand_dxydxy": (0, 0.001, 50, r"$\mathrm{cov}(d_{xy},d_{xy})$"),
    "pfcand_dzdz": (0, 0.005, 50, r"$\mathrm{cov}(d_z,d_z)$"),
    "pfcand_dphidphi": (0, 0.0001, 50, r"$\mathrm{cov}(\phi,\phi)$"),
    "pfcand_detadeta": (0, 0.0001, 50, r"$\mathrm{cov}(\eta,\eta)$"),
    "pfcand_dxydz": (-5e-4, 5e-4, 50, r"$\mathrm{cov}(d_{xy},d_z)$"),
    "pfcand_dphidxy": (-2e-4, 2e-4, 50, r"$\mathrm{cov}(\phi,d_{xy})$"),
    "pfcand_phidz": (-2e-4, 2e-4, 50, r"$\mathrm{cov}(\phi,d_z)$"),
    "pfcand_phictgtheta": (-1e-5, 1e-5, 50, r"$\mathrm{cov}(\phi,\cot\theta)$"),
    "pfcand_dxyctgtheta": (-1e-4, 1e-4, 50, r"$\mathrm{cov}(d_{xy},\cot\theta)$"),
    "pfcand_dlambdadz": (-1e-3, 1e-3, 50, r"$\mathrm{cov}(\lambda,d_z)$"),
    "pfcand_cctgtheta": (-1e-6, 1e-6, 50, r"$\mathrm{cov}(C,\cot\theta)$"),
    "pfcand_phic": (-1e-7, 1e-7, 50, r"$\mathrm{cov}(\phi,C)$"),
    "pfcand_dxyc": (-1e-5, 1e-5, 50, r"$\mathrm{cov}(d_{xy},C)$"),
    "pfcand_cdz": (-1e-5, 1e-5, 50, r"$\mathrm{cov}(C,d_z)$"),

    # -------------------------------------------------- b-tagging ----------
    "pfcand_btagSip2dVal": (-0.1, 0.5, 60, r"2D IP Value"),
    "pfcand_btagSip2dSig": (-30, 30, 60, r"2D IP Significance"),
    "pfcand_btagSip3dVal": (-0.1, 1.0, 60, r"3D IP Value"),
    "pfcand_btagSip3dSig": (-30, 30, 60, r"3D IP Significance"),
    "pfcand_btagJetDistVal": (0, 0.1, 50, r"Distance to Jet Axis"),
    "pfcand_btagJetDistSig": (-20, 20, 50, r"Distance to Jet Axis Significance"),

    # ------------------------------------------- dE/dx and PID (pads) ------
    "pfcand_dEdx_pads_type": (-10, 10, 20, r"dE/dx type (Pads)"),
    "pfcand_dEdx_pads_value": (-10, 20, 60, r"dE/dx (Pads)"),
    "pfcand_dEdx_pads_error": (-10, 5, 60, r"dE/dx error (Pads)"),
    "pfcand_PID_pval_pads_ele": (-1.2, 1.2, 48, r"Electron $p$-val (Pads)"),
    "pfcand_PID_pval_pads_mu": (-1.2, 1.2, 48, r"Muon $p$-val (Pads)"),
    "pfcand_PID_pval_pads_pi": (-1.2, 1.2, 48, r"Pion $p$-val (Pads)"),
    "pfcand_PID_pval_pads_kaon": (-1.2, 1.2, 48, r"Kaon $p$-val (Pads)"),
    "pfcand_PID_pval_pads_proton": (-1.2, 1.2, 48, r"Proton $p$-val (Pads)"),

    # ------------------------------------------ dE/dx and PID (wires) ------
    "pfcand_dEdx_wires_type": (-10, 10, 20, r"dE/dx type (Wires)"),
    "pfcand_dEdx_wires_value": (-10, 15, 50, r"dE/dx (Wires)"),
    "pfcand_dEdx_wires_error": (-10, 5, 60, r"dE/dx error (Wires)"),
    "pfcand_PID_pval_wires_ele": (-10, 1.5, 46, r"Electron $p$-val (Wires)"),
    "pfcand_PID_pval_wires_mu": (-10, 1.5, 46, r"Muon $p$-val (Wires)"),
    "pfcand_PID_pval_wires_pi": (-10, 1.5, 46, r"Pion $p$-val (Wires)"),
    "pfcand_PID_pval_wires_kaon": (-10, 1.5, 46, r"Kaon $p$-val (Wires)"),
    "pfcand_PID_pval_wires_proton": (-10, 1.5, 46, r"Proton $p$-val (Wires)"),

    # ------------------------------------------------------- thrust --------
    "EVT_Thrust_Mag": (0.5, 1.0, 50, r"Thrust $T$"),
    "EVT_Thrust_X": (-1, 1, 50, r"Thrust axis $x$"),
    "EVT_Thrust_Y": (-1, 1, 50, r"Thrust axis $y$"),
    "EVT_Thrust_Z": (-1, 1, 50, r"Thrust axis $z$"),
    "EVT_Thrust_cosTheta": (-1, 1, 50, r"$\cos\theta_{\rm thrust}$"),
    "EVT_Evis": (0, 120, 60, r"Visible energy $E_{\rm vis}$ [GeV]"),

    # ------------------------------------------------------- debug ---------
    "pfcand_dEdx_len": (0, 60, 60, r"dE/dx array length (jet 0)"),
    "pfcand_E_len": (0, 60, 60, r"Energy array length (jet 0)"),
    "pfcand_pval_ele_len": (0, 60, 60, r"$p$-val array length (jet 0)"),
}


# ---------------------------------------------------------------------------
# Finder modules (v0n_*, svn_*, phikk_*, dstar_*) - same tables as before.
# ---------------------------------------------------------------------------

def _add(prefix, table, label_prefix):
    for suffix, (lo, hi, nb, lab) in table.items():
        VAR_MAP[f"{prefix}{suffix}"] = (lo, hi, nb, f"{label_prefix} {lab}")


_LEG_COMMON = {
    "dEdx_pads_value": (-2, 20, 88, r"dE/dx (Pads)"),
    "dEdx_pads_error": (-2, 5, 70, r"dE/dx error (Pads)"),
    "dEdx_wires_value": (-2, 15, 68, r"dE/dx (Wires)"),
    "dEdx_wires_error": (-2, 5, 70, r"dE/dx error (Wires)"),
    "isChargedHad": (-1.5, 1.5, 3, r"PF charged-hadron label"),
}
_EXCL_TRK = {
    "q": (-1.5, 1.5, 3, r"charge"),
    "p": (0, 30, 60, r"$p$ [GeV]"),
    "costheta": (-1, 1, 50, r"$\cos\theta$"),
    "d0": (-0.3, 0.3, 60, r"$d_0$ [cm]"),
    "z0": (-2, 2, 50, r"$z_0$ [cm]"),
    "sigd0": (-20, 20, 60, r"$d_0/\sigma_{d_0}$"),
    "nvdet": (0, 5, 5, r"VDET hits"),
    "nitc": (0, 9, 9, r"ITC hits"),
    "chi2ndf": (0, 5, 50, r"track $\chi^2/$ndf"),
    "isprim": (0, 2, 2, r"primary-track flag"),
}

VAR_MAP["n_v0n_event"] = (0, 10, 10, r"$N_{V^0}$ (new) per event")
_add("v0n_", {
    "pdg": (0, 3500, 70, r"$|{\rm PDG}|$"),
    "invM": (0.3, 1.3, 100, r"invariant mass [GeV]"),
    "alpha": (-1, 1, 50, r"Armenteros $\alpha$"),
    "qt": (0, 0.3, 60, r"Armenteros $q_T$ [GeV]"),
    "chi2": (0, 20, 50, r"vertex $\chi^2$"),
    "dxyz": (0, 60, 60, r"$d_{3D}$ from PV [cm]"),
    "px": (-20, 20, 80, r"$p_x$ [GeV]"),
    "py": (-20, 20, 80, r"$p_y$ [GeV]"),
    "pz": (-40, 40, 80, r"$p_z$ [GeV]"),
    "cosPointing": (0.9, 1.0, 50, r"$\cos\alpha_{\rm pointing}$"),
    "pointSig": (0, 10, 50, r"pointing significance"),
    "tight": (0, 2, 2, r"tight tier"),
    "bandSig": (-6, 6, 60, r"band pull"),
    "massSig": (-6, 6, 60, r"mass pull"),
    "vx": (-40, 40, 80, r"vertex $x$ [cm]"),
    "vy": (-40, 40, 80, r"vertex $y$ [cm]"),
    "vz": (-80, 80, 80, r"vertex $z$ [cm]"),
    "cov_xx": (0, 0.05, 50, r"cov$(x,x)$ [cm$^2$]"),
    "cov_yy": (0, 0.05, 50, r"cov$(y,y)$ [cm$^2$]"),
    "cov_zz": (0, 0.05, 50, r"cov$(z,z)$ [cm$^2$]"),
    "jetIdx": (-1.5, 1.5, 3, r"jet index"),
    "z": (0, 1, 50, r"$z = p/p_{\rm jet}$"),
    "zL": (0, 1, 50, r"$z_L$"),
    "ptRel": (0, 2, 50, r"$p_{T,\rm rel}$ [GeV]"),
    "dRjet": (0, 1, 50, r"$\Delta R$ to jet"),
    "rankInJet": (0, 6, 6, r"rank in jet"),
    "Lxy": (0, 40, 80, r"$L_{xy}$ [cm]"),
    "LxySig": (0, 200, 80, r"$L_{xy}/\sigma$"),
    "LxyzSig": (0, 200, 80, r"$L_{xyz}/\sigma$"),
    "baryon": (-1.5, 1.5, 3, r"baryon number"),
    "nShared": (0, 3, 3, r"shared daughters"),
    "svnCosPoint": (-1, 1, 50, r"$\cos$ pointing to SV"),
    "svnPointSig": (0, 20, 50, r"pointing significance to SV"),
}, r"$V^0$")
for _t, _tl in (("trk1", "leg 1"), ("trk2", "leg 2")):
    _add(f"v0n_{_t}_", dict(_LEG_COMMON, **{
        "q": (-1.5, 1.5, 3, r"charge"),
        "p": (0, 30, 60, r"$p$ [GeV]"),
        "nTPC": (0, 25, 25, r"TPC hits"),
    }), rf"$V^0$ {_tl}")

VAR_MAP["n_svn_event"] = (0, 10, 10, r"$N_{\rm SV}$ (new) per event")
_add("svn_", {
    "mass": (0, 6, 60, r"mass [GeV]"),
    "chi2": (0, 20, 50, r"vertex $\chi^2$"),
    "dxyz": (0, 5, 50, r"$d_{3D}$ from PV [cm]"),
    "dx": (-3, 3, 60, r"$\Delta x$ from PV [cm]"),
    "dy": (-3, 3, 60, r"$\Delta y$ from PV [cm]"),
    "dz": (-3, 3, 60, r"$\Delta z$ from PV [cm]"),
    "p": (0, 40, 60, r"$p$ [GeV]"),
    "cosPointing": (0.9, 1.0, 50, r"$\cos\alpha_{\rm pointing}$"),
    "pointSig": (0, 20, 50, r"pointing significance"),
    "ntracks": (2, 10, 8, r"track multiplicity"),
    "sigL": (0, 50, 50, r"$L/\sigma_L$"),
    "cov_xx": (0, 0.01, 50, r"cov$(x,x)$ [cm$^2$]"),
    "cov_yy": (0, 0.01, 50, r"cov$(y,y)$ [cm$^2$]"),
    "cov_zz": (0, 0.01, 50, r"cov$(z,z)$ [cm$^2$]"),
}, r"SV$_{\rm new}$")

VAR_MAP["n_phikk_event"] = (0, 10, 10, r"$N_{\phi\to KK}$ per event")
_add("phikk_", {
    "invM": (0.98, 1.10, 60, r"$m_{KK}$ [GeV]"),
    "p": (0, 40, 60, r"$p$ [GeV]"),
    "px": (-20, 20, 80, r"$p_x$ [GeV]"),
    "py": (-20, 20, 80, r"$p_y$ [GeV]"),
    "pz": (-40, 40, 80, r"$p_z$ [GeV]"),
    "alpha": (-1, 1, 50, r"Armenteros $\alpha$"),
    "qt": (0, 0.15, 60, r"Armenteros $q_T$ [GeV]"),
    "bandEll": (0, 5, 50, r"band ellipse"),
    "chi2": (0, 20, 50, r"vertex $\chi^2$"),
    "vx": (-1, 1, 50, r"vertex $x$ [cm]"),
    "vy": (-1, 1, 50, r"vertex $y$ [cm]"),
    "vz": (-5, 5, 50, r"vertex $z$ [cm]"),
    "dpv": (0, 1, 50, r"distance to PV [cm]"),
    "dpvSig": (0, 20, 50, r"distance to PV / $\sigma$"),
    "same_sign": (0, 2, 2, r"same-sign pair"),
    "wp": (0, 4, 4, r"working point"),
    "tight": (0, 2, 2, r"tight"),
}, r"$\phi\to KK$")
for _t, _tl in (("trk1", "leg 1"), ("trk2", "leg 2")):
    _add(f"phikk_{_t}_", dict(_LEG_COMMON, **_EXCL_TRK), rf"$\phi\to KK$ {_tl}")

VAR_MAP["n_dstar_event"] = (0, 10, 10, r"$N_{D^*}$ per event")
VAR_MAP["n_d0fits_event"] = (0, 500, 50, r"$D^0$ two-track fits per event")
_add("dstar_", {
    "m_kpi": (1.7, 2.0, 60, r"$m_{K\pi}$ [GeV]"),
    "dm": (0.139, 0.165, 52, r"$\Delta m = m_{K\pi\pi} - m_{K\pi}$ [GeV]"),
    "p": (0, 50, 50, r"$p$ [GeV]"),
    "px": (-30, 30, 60, r"$p_x$ [GeV]"),
    "py": (-30, 30, 60, r"$p_y$ [GeV]"),
    "pz": (-45, 45, 60, r"$p_z$ [GeV]"),
    "costheta": (-1, 1, 50, r"$\cos\theta$"),
    "xE": (0, 1, 50, r"$x_E$"),
    "chi2": (0, 20, 50, r"$D^0$ vertex $\chi^2$"),
    "vx": (-1, 1, 50, r"$D^0$ vertex $x$ [cm]"),
    "vy": (-1, 1, 50, r"$D^0$ vertex $y$ [cm]"),
    "vz": (-5, 5, 50, r"$D^0$ vertex $z$ [cm]"),
    "dpv": (0, 1, 50, r"$D^0$ distance to PV [cm]"),
    "dpvSig": (0, 30, 60, r"$D^0$ distance to PV / $\sigma$"),
    "cosPoint": (-1, 1, 50, r"$D^0$ $\cos$ pointing"),
    "cosThetaStar": (-1, 1, 50, r"$\cos\theta^*$"),
    "rs": (0, 2, 2, r"right-sign"),
    "loose": (0, 2, 2, r"loose"),
    "tight": (0, 2, 2, r"tight"),
    "nsec": (0, 4, 4, r"secondary legs"),
}, r"$D^*$")
for _pfx, _tl in (("dstar_trkK_", r"$D^*$ $K$"), ("dstar_trkPi_", r"$D^*$ $\pi$"),
                  ("dstar_trkPis_", r"$D^*$ $\pi_{\rm slow}$")):
    _add(_pfx, dict(_LEG_COMMON, **_EXCL_TRK,
                    pool=(0, 2, 2, r"pool class (0 prim / 1 sec)")), _tl)


# ======================================================================
# 5. WORKER-SIDE CONFIGURATION
# ======================================================================

_BINNING = {}
_CFG = {}


def _init_worker(binning, cfg):
    global _BINNING, _CFG
    _BINNING = binning
    _CFG = cfg


_TRANSFORMS = {
    None: lambda x: x,
    "abs": lambda x: abs(x),
    "cos": lambda x: np.cos(x),
    "abscos": lambda x: abs(np.cos(x)),
    "tanh": lambda x: np.tanh(x),
    "abstanh": lambda x: abs(np.tanh(x)),
    "log": lambda x: np.log(x),
    "log10": lambda x: np.log10(x),
}


def preselection_branches(presel, keys):
    return [b for level in presel for b in presel[level] if b in keys]


def _bounds_mask(values, spec):
    mask = None
    if "min" in spec:
        mask = values > spec["min"]
    if "max" in spec:
        upper = values < spec["max"]
        mask = upper if mask is None else mask & upper
    if "equals" in spec:
        eq = values == spec["equals"]
        mask = eq if mask is None else mask & eq
    return mask


def _cut_values(chunk, branch, spec, keys):
    if branch not in keys:
        return None
    return _TRANSFORMS[spec.get("transform")](chunk[branch])


def _event_mask(chunk, keys, cfg):
    presel = cfg["presel"]
    mask = None

    for branch, spec in presel.get("event", {}).items():
        values = _cut_values(chunk, branch, spec, keys)
        if values is None:
            continue
        if values.ndim > 1:
            values = ak.firsts(values)
        sub = _bounds_mask(values, spec)
        if sub is not None:
            mask = sub if mask is None else mask & sub

    for branch, spec in presel.get("jet", {}).items():
        values = _cut_values(chunk, branch, spec, keys)
        if values is None:
            continue
        sub = _bounds_mask(values, spec)
        if sub is None:
            continue
        reduce_fn = ak.any if spec.get("require") == "any" else ak.all
        sub = reduce_fn(sub, axis=1) & (ak.num(values, axis=1) > 0)
        mask = sub if mask is None else mask & sub

    if mask is None:
        return np.ones(len(chunk), dtype=bool)
    return ak.fill_none(mask, False)


def _object_masks(chunk, keys, cfg):
    masks = {}
    for level, prefix in OBJECT_PREFIXES.items():
        combined = None
        for branch, spec in cfg["presel"].get(level, {}).items():
            values = _cut_values(chunk, branch, spec, keys)
            if values is None:
                continue
            sub = _bounds_mask(values, spec)
            if sub is not None:
                combined = sub if combined is None else combined & sub
        if combined is not None:
            masks[prefix] = combined
    return masks


def _extract_values(chunk, var, obj_masks=None):
    arr = chunk[var]
    for prefix, mask in (obj_masks or {}).items():
        if var.startswith(prefix) and arr.ndim == mask.ndim:
            arr = arr[mask]
            break
    vals = ak.flatten(arr, axis=None) if arr.ndim > 1 else arr
    if len(vals) == 0:
        return np.empty(0)
    try:
        return np.asarray(vals)
    except Exception:
        return np.asarray(ak.to_numpy(ak.fill_none(vals, np.nan)), dtype=float)


def _prepare_tree(tree, var_names, cfg):
    keys = set(tree.keys())
    present = [v for v in var_names if v in keys]
    if not present:
        return None, None, keys
    extra = preselection_branches(cfg["presel"], keys)
    return present, list(dict.fromkeys(present + extra)), keys


def _hist_worker(args):
    file_path, var_names = args
    counts = {v: np.zeros(_BINNING[v][2], dtype=np.float64) for v in var_names}
    n_selected = 0

    try:
        with uproot.open(file_path) as root_file:
            if "events" not in root_file:
                return counts, 0
            tree = root_file["events"]
            present, branches, keys = _prepare_tree(tree, var_names, _CFG)
            if not present:
                return counts, 0

            for chunk in tree.iterate(branches, step_size=_CFG["step"], library="ak"):
                chunk = chunk[_event_mask(chunk, keys, _CFG)]
                if len(chunk) == 0:
                    continue
                obj_masks = _object_masks(chunk, keys, _CFG)
                n_selected += len(chunk)
                for var in present:
                    vals = _extract_values(chunk, var, obj_masks)
                    if vals.size == 0:
                        continue
                    xmin, xmax, nbins, _ = _BINNING[var]
                    hist, _ = np.histogram(vals, bins=nbins, range=(xmin, xmax))
                    counts[var] += hist

    except Exception as exc:
        print(f"  [warn] {os.path.basename(file_path)}: {exc}")

    return counts, n_selected


def _probe_worker(args):
    file_path, var_names = args
    sample = {}

    try:
        with uproot.open(file_path) as root_file:
            if "events" not in root_file:
                return sample
            tree = root_file["events"]
            present, branches, keys = _prepare_tree(tree, var_names, _CFG)
            if not present:
                return sample

            rng = np.random.default_rng(0)
            for chunk in tree.iterate(branches, step_size="50 MB", library="ak"):
                chunk = chunk[_event_mask(chunk, keys, _CFG)]
                if len(chunk) == 0:
                    break
                obj_masks = _object_masks(chunk, keys, _CFG)
                for var in present:
                    vals = _extract_values(chunk, var, obj_masks)
                    if vals.size > 20000:
                        vals = vals[rng.choice(vals.size, 20000, replace=False)]
                    if vals.size:
                        sample[var] = vals
                break

    except Exception:
        pass

    return sample


# ======================================================================
# 6. FILE DISCOVERY
# ======================================================================

def get_root_files(directory, max_files=None):
    files = sorted(glob.glob(os.path.join(directory, "**", "*.root"), recursive=True))
    if not files:
        print(f"WARNING: no ROOT files found in {directory}")
    if max_files:
        files = files[:max_files]
    return files


def get_flavor_files(prod_dir, flavor, max_files=None):
    """<prod_dir>/<flavor>/**/*.root (chunked) or <prod_dir>/<flavor>.root."""
    directory = os.path.join(prod_dir, flavor)
    if os.path.isdir(directory):
        return get_root_files(directory, max_files)
    single = directory + ".root"
    if os.path.isfile(single):
        return [single]
    print(f"WARNING: no ROOT files for flavour {flavor} under {prod_dir}")
    return []


# ======================================================================
# 7. FILL ONE SAMPLE
# ======================================================================

def fill_one_sample(files, var_list, executor, tag):
    """Accumulate {var: counts} over every file of one (production, flavour)."""
    totals = {v: np.zeros(VAR_MAP[v][2], dtype=np.float64) for v in var_list}
    if not files:
        print(f"    {tag:20s}: no files")
        return totals, 0

    tasks = [(p, tuple(var_list)) for p in files]
    chunksize = max(1, len(files) // (NUM_WORKERS * 4))
    n_selected = 0
    t0 = time.time()

    for counts, nsel in executor.map(_hist_worker, tasks, chunksize=chunksize):
        for var, arr in counts.items():
            totals[var] += arr
        n_selected += nsel

    elapsed = time.time() - t0
    rate = len(files) / elapsed if elapsed > 0 else 0.0
    print(f"    {tag:20s}: {len(files):4d} files, {n_selected:>12,} selected events "
          f"({elapsed:6.1f} s, {rate:5.1f} files/s)")
    return totals, n_selected


# ======================================================================
# 8. PLOTTING
# ======================================================================

def default_cut_label(branch, spec):
    name = branch.replace("_", r"\_")
    name = (r"\texttt{%s}" % name) if HAS_LATEX else branch
    if spec.get("transform") in ("abs", "abscos", "abstanh"):
        name = "|%s|" % name
    parts = []
    if "min" in spec:
        parts.append("$> %g$" % spec["min"])
    if "max" in spec:
        parts.append("$< %g$" % spec["max"])
    if "equals" in spec:
        parts.append("$= %g$" % spec["equals"])
    return f"{name} " + " ".join(parts)


def cut_labels(max_lines=4):
    labels = [spec.get("label") or default_cut_label(branch, spec)
              for level in ("event", "jet", "constituent", "sv", "v0")
              for branch, spec in PRESELECTION.get(level, {}).items()]
    if len(labels) > max_lines:
        labels = labels[:max_lines - 1] + ["+ %d more cuts" % (len(labels) - max_lines + 1)]
    return labels


def _normalize_to_unit_area(raw):
    """Return raw counts scaled so that sum(counts) * binwidth == 1."""
    total = raw.sum()
    if total <= 0:
        return None, 0.0
    return raw / total, float(total)


def draw_flavour_variable(flavour, var, counts_by_prod, scale, outdir):
    xmin, xmax, nbins, label = VAR_MAP[var]
    edges = np.linspace(xmin, xmax, nbins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    binwidth = edges[1] - edges[0]

    series = []          # [(prod_key, prod_label, colour, normalized, n_entries)]
    for prod_key, prod_label, _, colour in PRODUCTIONS:
        raw = counts_by_prod.get(prod_key, {}).get(var)
        if raw is None or raw.sum() == 0:
            continue
        norm, n_entries = _normalize_to_unit_area(raw)
        if norm is None:
            continue
        series.append((prod_key, prod_label, colour, norm, n_entries))

    if len(series) < 1:
        return False

    ref_key = series[0][0]
    non_ref = [s for s in series[1:] if s[0] != ref_key]
    n_ratio = len(non_ref)

    if n_ratio == 0:
        fig, ax1 = plt.subplots(figsize=(6, 4.5))
        axes = [ax1]
    else:
        fig, axes = plt.subplots(
            1 + n_ratio, 1, sharex=True,
            figsize=(6, 4 + 1.2 * n_ratio),
            gridspec_kw={"height_ratios": [3] + [1] * n_ratio},
        )
    ax1 = axes[0]
    plt.subplots_adjust(hspace=0.07, left=0.15, right=0.97,
                        top=0.90, bottom=0.12)

    header = r"\textbf{ALEPH} MC comparison" if HAS_LATEX else "ALEPH MC comparison"
    ax1.text(0.0, 1.02, header, transform=ax1.transAxes, ha="left", va="bottom")
    ax1.text(1.0, 1.02, FLAVOR_LABELS.get(flavour, flavour),
             transform=ax1.transAxes, ha="right", va="bottom")

    for i, text in enumerate(cut_labels()):
        ax1.text(0.03, 0.95 - 0.07 * i, text, transform=ax1.transAxes,
                 ha="left", va="top", fontsize=8)

    ref_norm = next(s[3] for s in series if s[0] == ref_key)

    for prod_key, prod_label, colour, norm, n_entries in series:
        ax1.stairs(norm, edges, color=colour, fill=True, alpha=0.13, linewidth=0)
        ax1.stairs(norm, edges, color=colour, linewidth=1.6,
                   label=f"{prod_label}  ($N = {int(n_entries):,}$)")

    ax1.set_ylabel("Normalised to unit area")
    ax1.legend(fontsize=8, loc="upper right", frameon=False)

    # ---- y-limits
    ymax = max((s[3].max() for s in series), default=1.0)
    if ymax <= 0:
        ymax = 1.0

    if scale == "log":
        ax1.set_yscale("log")
        # smallest positive value across all curves
        positives = [v for s in series for v in s[3] if v > 0]
        ymin = min(positives) * 0.3 if positives else 1e-6
        ax1.set_ylim(ymin, ymax * 5)
    else:
        ax1.set_yscale("linear")
        ax1.set_ylim(0, ymax * 1.35)

    # ---- ratio panels (each non-reference production / reference)
    for i, (prod_key, prod_label, colour, norm, _) in enumerate(non_ref):
        ax = axes[i + 1]
        safe = ref_norm > 0
        ratio = np.zeros(nbins)
        ratio[safe] = norm[safe] / ref_norm[safe]
        ax.step(centers, ratio, where="mid", color=colour, lw=1.2)
        ax.axhline(1.0, color="red", lw=0.8, ls="--")
        ax.set_ylabel(f"{prod_key}/{ref_key}", fontsize=8)
        ax.set_ylim(0.5, 1.5)
        ax.grid(axis="y", alpha=0.2)
        if i == len(non_ref) - 1:
            ax.set_xlabel(label)

    if n_ratio == 0:
        ax1.set_xlabel(label)

    ax1.set_xlim(xmin, xmax)
    os.makedirs(outdir, exist_ok=True)
    fig.savefig(os.path.join(outdir, f"{var}.png"), bbox_inches="tight", dpi=DPI)
    plt.close(fig)
    return True


# ======================================================================
# 9. AUTO-RANGE (optional)
# ======================================================================

def _auto_label(var):
    if var in VAR_MAP and VAR_MAP[var][3]:
        return VAR_MAP[var][3]
    escaped = var.replace("_", r"\_")
    return (r"\texttt{%s}" % escaped) if HAS_LATEX else var


def _write_range_file(binning, var_list):
    try:
        with open(RANGE_CACHE, "w") as handle:
            handle.write("# Auto-derived binning - paste into VAR_MAP if it looks right.\n")
            handle.write("VAR_MAP = {\n")
            for var in var_list:
                lo, hi, nb, label = binning[var]
                handle.write(f'    "{var}": ({lo:.6g}, {hi:.6g}, {nb}, r"{label}"),\n')
            handle.write("}\n")
        print(f"Suggested binning written to {RANGE_CACHE}")
    except Exception as exc:
        print(f"Could not write {RANGE_CACHE}: {exc}")


def probe_ranges(var_list, all_files, executor, n_probe):
    print("\nSTAGE 3: RANGE PROBE (pooled across all productions and flavours)")
    tasks = []
    for flavour in FLAVORS:
        for prod_key, _, _, _ in PRODUCTIONS:
            for path in all_files.get(flavour, {}).get(prod_key, [])[:n_probe]:
                tasks.append((path, tuple(var_list)))

    pooled = {}
    for sample in executor.map(_probe_worker, tasks):
        for var, vals in sample.items():
            pooled.setdefault(var, []).append(vals)

    binning = dict(VAR_MAP)
    for var in var_list:
        chunks = pooled.get(var)
        if not chunks:
            continue
        vals = np.concatenate(chunks)
        vals = vals[np.isfinite(vals)]
        if vals.size < 100:
            continue

        lo = float(np.percentile(vals, 0.2))
        hi = float(np.percentile(vals, 99.8))
        if hi <= lo:
            lo, hi = float(vals.min()), float(vals.max())
        if hi <= lo:
            hi = lo + 1.0

        if np.allclose(vals, np.round(vals)) and (hi - lo) <= 100:
            lo, hi = float(np.floor(lo)), float(np.ceil(hi) + 1)
            nbins = int(hi - lo)
        else:
            pad = 0.05 * (hi - lo)
            lo, hi = lo - pad, hi + pad
            nbins = 50

        binning[var] = (lo, hi, max(nbins, 2), _auto_label(var))

    _write_range_file(binning, var_list)
    return binning


# ======================================================================
# 10. MAIN
# ======================================================================

def report_preselection(all_files):
    print("\nPRESELECTION")
    probes = {}
    for flavour in FLAVORS:
        for prod_key, _, _, _ in PRODUCTIONS:
            files = all_files.get(flavour, {}).get(prod_key, [])
            if files:
                probes[f"{flavour}/{prod_key}"] = files[0]

    keys = {}
    for tag, path in probes.items():
        try:
            with uproot.open(path) as f:
                keys[tag] = set(f["events"].keys()) if "events" in f else set()
        except Exception as exc:
            print(f"  could not read {tag} probe file: {exc}")
            keys[tag] = set()

    # check just the first flavour/production pair for brevity
    sample_tags = list(keys.keys())[:3]
    for level in ("event", "jet", "constituent", "sv", "v0"):
        cuts = PRESELECTION.get(level, {})
        if not cuts:
            continue
        scope = ("removes events" if level in ("event", "jet")
                 else f"filters {OBJECT_PREFIXES[level]}* objects")
        print(f"  [{level}] ({scope})")
        for branch, spec in cuts.items():
            bits = [f"{tag}: {'ok' if branch in keys[tag] else 'MISSING'}"
                    for tag in sample_tags]
            label = spec.get("label") or default_cut_label(branch, spec)
            print(f"      {branch:24s} {'  '.join(bits):40s} {label}")
    print()


def parse_args():
    p = argparse.ArgumentParser(
        description="ALEPH MC comparison: Paper vs HVFL vs Pythia 8.3 per flavour")
    p.add_argument("--workers", type=int, default=NUM_WORKERS)
    p.add_argument("--auto-range", action="store_true",
                   help="derive histogram ranges from the files themselves")
    p.add_argument("--probe-files", type=int, default=3,
                   help="files per (flavour, production) used by --auto-range")
    p.add_argument("--outdir", default=None, help="base output directory")
    p.add_argument("--vars", nargs="+", default=None, help="subset of variables")
    p.add_argument("--flavors", nargs="+", default=None,
                   help="subset of flavours (default: Zbb Zcc Zss Zuu Zdd)")
    p.add_argument("--max-files", type=int, default=None,
                   help="cap files per (flavour, production) for quick tests")
    p.add_argument("--step-size", default=STEP_SIZE,
                   help="uproot read step, e.g. '50 MB' if memory is tight")
    return p.parse_args()


def main():
    global NUM_WORKERS, OUTPUT_BASE

    args = parse_args()
    NUM_WORKERS = args.workers
    OUTPUT_BASE = args.outdir or OUTPUT_BASE

    var_list = args.vars if args.vars else list(VAR_MAP)
    for v in [v for v in var_list if v not in VAR_MAP]:
        print(f"Skipping unknown variable: {v}")
    var_list = [v for v in var_list if v in VAR_MAP]

    flavour_list = args.flavors if args.flavors else FLAVORS
    flavour_list = [f for f in FLAVORS if f in flavour_list]  # preserve order

    cfg = {"presel": PRESELECTION, "step": args.step_size}

    start = time.time()
    print("=" * 70)
    print(f"ALEPH MC COMPARISON (workers: {NUM_WORKERS})")
    print(f"Productions: {', '.join(p[0] for p in PRODUCTIONS)}")
    print(f"Flavours   : {', '.join(flavour_list)}")
    print(f"Variables  : {len(var_list)}")
    print(f"Output     : {OUTPUT_BASE}/<flavour>/{{linear,log}}/")
    print("=" * 70)

    # ---- Stage 0: file discovery  (grouped per flavour first)
    all_files = {}
    for flavour in flavour_list:
        all_files[flavour] = {}
        for prod_key, _, prod_dir, _ in PRODUCTIONS:
            all_files[flavour][prod_key] = get_flavor_files(
                prod_dir, flavour, args.max_files)
        counts = ", ".join(f"{p[0]}:{len(all_files[flavour][p[0]])}"
                           for p in PRODUCTIONS)
        print(f"  {flavour}: {counts}")

    report_preselection(all_files)

    # ---- optional: auto-range
    if args.auto_range:
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=NUM_WORKERS,
            initializer=_init_worker,
            initargs=(VAR_MAP, cfg),
        ) as executor:
            binning = probe_ranges(var_list, all_files, executor, args.probe_files)
    else:
        binning = dict(VAR_MAP)

    # ---- Stage 1-3: per flavour, fill the three productions, then draw
    total_made = 0
    with concurrent.futures.ProcessPoolExecutor(
        max_workers=NUM_WORKERS,
        initializer=_init_worker,
        initargs=(binning, cfg),
    ) as executor:
        for flavour in flavour_list:
            print("\n" + "=" * 70)
            print(f"FLAVOUR {flavour}  ({FLAVOR_LABELS.get(flavour, flavour)})")
            print("=" * 70)

            counts_by_prod = {}
            for prod_key, prod_label, _, _ in PRODUCTIONS:
                files = all_files[flavour][prod_key]
                counts, _ = fill_one_sample(
                    files, var_list, executor, f"{prod_key}")
                counts_by_prod[prod_key] = counts

            print(f"  drawing {len(var_list)} variables x {len(SCALES)} scales ...")
            made = 0
            for var in var_list:
                for scale in SCALES:
                    outdir = os.path.join(OUTPUT_BASE, flavour, scale)
                    try:
                        if draw_flavour_variable(
                                flavour, var, counts_by_prod, scale, outdir):
                            made += 1
                    except Exception as exc:
                        print(f"    ERROR drawing {var} ({scale}): {exc}")
            total_made += made
            print(f"  -> {made} plots written under {os.path.join(OUTPUT_BASE, flavour)}")

    print("\n" + "=" * 70)
    print(f"TOTAL PLOTS: {total_made}")
    print(f"TOTAL TIME : {time.time() - start:.1f} s")
    print("=" * 70)


if __name__ == "__main__":
    main()
