#!/usr/bin/env python3
"""Synthetic EDM4hep input for the stage1 CI (no ALEPH data involved).

Writes one data-like and one MC-like podio ROOT file with the collection names and
types of the converted ALEPH files that stage1.py reads:

  <out>/1994/synthetic_data.root   no MC truth (like the data)
  <out>/QQB/synthetic_mc.root      plus MCParticles, trackMCLink, ecalClusterMCLink, GeneratorInfo
  <out>/beamspot.json              per-run beam-spot file for $ALEPH_BEAMSPOT_JSON (made-up values)

The events are random but deterministic: event i of a file depends only on (--seed, data or
mc, i), and the files are the same for any number of worker processes (--jobs; the events
are written in chunks of CHUNK events and joined in order). Z -> q qbar at the Z pole with
prompt charged hadrons, photons, neutral hadrons, leptons, K0s -> pi pi, Lambda -> p pi,
phi -> K K, D*+ -> D0 pi (D0 -> K pi), D+ -> K pi pi and B -> D (-> K pi pi) + 3 pi with
displaced vertices, measured as perigee helices with covariances, hit counts and dE/dx.
The data events are spread over three selected runs of the run list, each with its own
beam-spot centre, and one rejected run; a few events fail the hadronic class bit.

The input is meant to exercise the code paths of stage1, not to look like ALEPH data:
the detector geometry is the public one and the track resolutions are of the published
size, kinematics and fractions are generic; nothing is fitted to ALEPH data.
  - Track quality: some tracks have degraded hit efficiencies or a short TPC segment
    (broad VDET/ITC/TPC hit counts, some below the TPC-hit cut), d0/z0 residuals with
    non-Gaussian tails that the covariance does not describe, a chi2 tail across the
    chi2/ndf cut, a covariance that is not positive definite, or |z0| around the cut;
    some particles have a second, short track segment without an energy-flow object,
    which mostly fails the track selection.
  - Awkward events (AWKWARD_FRACTION of the hadronic events; legal, but unusual): one
    track or none for the primary-vertex fit (pv_trivial = 1); a track with vanishing
    errors (singular vertex fit: pv_converged = 0); three tracks whose fit converges but
    whose refit after pruning is singular (pv_split_converged = 0 with pv_converged = 1);
    a nearly identical duplicate of a track (a collinear pair); a track with a huge error.
    The vertex-fit outcomes follow from the fit's eigenvalue floor in src/analyzer_pvnew.h.

Conventions (those of the converted files):
  lengths cm, momenta GeV, Bz = 1.5 T; one track state per track, perigee to the origin;
  stored D0/omega in the ALEPH sign convention (stage1 flips both): the point of
  closest approach is D0*(sin phi, -cos phi) and dphi/ds = omega = -q*0.0029979*Bz/pT;
  covariance packed lower-triangular in (d0, phi, omega, z0, tanLambda, time);
  subdetectorHitNumbers = (VDET, ITC, TPC); ParticleID index-parallel to RecoParticles,
  with the energy-flow type (0 charged hadron, 1 electron, 2 muon, 4 photon, 5 neutral hadron);
  a failed dE/dx leg carries the track omega as its value.
"""
import argparse
import json
import math
import multiprocessing
import os
import random
import shutil
import sys
import time

import numpy as np

KB = 0.0029979 * 1.5          # pT [GeV] = KB / |omega [1/cm]|
E_BEAM = 45.6
CLASS_HADRONIC_BIT = 15       # class 16 = bit 15
CLASS_LEPTONPAIR_BIT = 14     # class 15 = bit 14
CHUNK = 200                   # events per worker task: fixed, so that the files do not depend on --jobs
DATA_YEAR = "1994"            # src/stage1.py reads the data from <input>/1994/ (its data process "1994")
REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")   # for src/run_list.py

# data runs: three selected runs and one rejected run of the run list, with the fraction of
# the data events and a made-up beam-spot centre [cm] (away from the origin) for each
DATA_RUN_FRACTIONS = (0.45, 0.30, 0.15, 0.10)
DATA_BEAMSPOTS = ((0.050, 0.020, 0.10), (0.060, 0.012, -0.25), (0.037, 0.027, 0.40), (0.055, 0.018, 0.0))

# event mix (generic, not tuned)
P_MUMU = 0.04                 # mu+ mu- events instead of hadronic ones (they fail the hadronic class bit)
FLAVOUR_CUMULATIVE = (0.15, 0.27, 0.42, 0.62)   # hadronic events: d 15 %, u 12 %, s 15 %, c 20 %, the rest b

# track quality, as fractions of the reconstructed tracks (generic, not tuned)
P_SHORT = 0.06        # short TPC segment (the track is lost early in the TPC): few TPC hits
P_DEGRADED = 0.30     # hit efficiencies drawn per track: broad hit-count distributions
P_TAIL = 0.08         # d0/z0 residuals scaled by 1 + Exp(mean 3), not reflected in the covariance
P_CHI2_TAIL = 0.10    # chi2 scaled by a log-uniform factor in [1, 20]
P_BAD_COV = 0.01      # covariance not positive definite: fails the track selection
P_FAR_Z0 = 0.01       # |z0| between 40 and 80 cm, around the |z0| <= 50 cm cut
P_SEGMENT = 0.25      # an extra short segment of the same particle, without an energy-flow object
P_NO_EF = 0.05        # a track without an energy-flow object (as in the converted files)
P_DEDX_FAIL = (0.15, 0.05)           # failed dE/dx measurement (pads, wires) of a track
P_DEDX_FAIL_SEGMENT = (0.6, 0.5)     # the same, for an extra short segment
AWKWARD_FRACTION = 0.025
AWKWARD_KINDS = ("one_track", "no_track", "precise", "split_singular", "duplicate", "huge_error")

M = {"pi": 0.13957039, "K": 0.493677, "p": 0.93827208, "e": 0.000511, "mu": 0.1056584,
     "K0S": 0.497611, "Lambda": 1.115683, "phi": 1.019461, "D0": 1.86484, "Dstar": 2.01026,
     "Dplus": 1.86966, "B": 5.27966, "Kstar0": 0.89555, "rho": 0.77526, "a1": 1.23,
     "gamma": 0.0, "nh": 0.497611}
PDG = {"pi": 211, "K": 321, "p": 2212, "e": 11, "mu": 13, "K0S": 310, "Lambda": 3122,
       "phi": 333, "D0": 421, "Dstar": 413, "Dplus": 411, "B": 511, "Kstar0": 313,
       "rho": 113, "a1": 20213, "gamma": 22, "nh": 130}
CTAU = {"K0S": 2.6844, "Lambda": 7.89, "D0": 0.01229, "Dplus": 0.03118, "B": 0.0455}  # cm

# helix-fit correlations in (d0, phi, omega, z0, tanLambda); positive definite
CORR = np.array([[1.0, 0.6, 0.3, 0.0, 0.0],
                 [0.6, 1.0, 0.5, 0.0, 0.0],
                 [0.3, 0.5, 1.0, 0.0, 0.0],
                 [0.0, 0.0, 0.0, 1.0, -0.6],
                 [0.0, 0.0, 0.0, -0.6, 1.0]])
VDET_LAYERS = ((6.3, 0.85), (10.8, 0.69))       # (radius cm, |cos theta| acceptance)
ITC_RADII = tuple(16.0 + 10.0 * k / 7 for k in range(8))
TPC_RADII = tuple(31.0 + 149.0 * k / 20 for k in range(21))
TPC_HALF_Z = 220.0
CLEAN = (1.0, 1.0, 1.0, 1e9)  # hit efficiencies (VDET, ITC, TPC) and TPC end radius of a clean track


# ----------------------------------------------------------------------------------------
# kinematics
# ----------------------------------------------------------------------------------------
def unit(v):
    n = math.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)
    return (v[0] / n, v[1] / n, v[2] / n)


def boost(p4, beta):
    """Lorentz boost of p4 = (px, py, pz, E) by velocity beta (3-vector)."""
    b2 = beta[0] ** 2 + beta[1] ** 2 + beta[2] ** 2
    if b2 <= 0.0:
        return p4
    gamma = 1.0 / math.sqrt(1.0 - b2)
    bp = beta[0] * p4[0] + beta[1] * p4[1] + beta[2] * p4[2]
    g2 = (gamma - 1.0) / b2
    return (p4[0] + g2 * bp * beta[0] + gamma * beta[0] * p4[3],
            p4[1] + g2 * bp * beta[1] + gamma * beta[1] * p4[3],
            p4[2] + g2 * bp * beta[2] + gamma * beta[2] * p4[3],
            gamma * (p4[3] + bp))


def two_body(rng, p4, m1, m2):
    """Isotropic two-body decay of p4 (mass from p4) into masses m1, m2."""
    mm = math.sqrt(max(p4[3] ** 2 - p4[0] ** 2 - p4[1] ** 2 - p4[2] ** 2, 0.0))
    q = math.sqrt(max((mm ** 2 - (m1 + m2) ** 2) * (mm ** 2 - (m1 - m2) ** 2), 0.0)) / (2 * mm)
    c = rng.uniform(-1.0, 1.0)
    s = math.sqrt(1.0 - c * c)
    ph = rng.uniform(0.0, 2 * math.pi)
    d = (q * s * math.cos(ph), q * s * math.sin(ph), q * c)
    beta = (p4[0] / p4[3], p4[1] / p4[3], p4[2] / p4[3])
    a = boost((d[0], d[1], d[2], math.sqrt(q * q + m1 * m1)), beta)
    b = boost((-d[0], -d[1], -d[2], math.sqrt(q * q + m2 * m2)), beta)
    return a, b


def p4_of(mass, p3):
    return (p3[0], p3[1], p3[2], math.sqrt(mass ** 2 + p3[0] ** 2 + p3[1] ** 2 + p3[2] ** 2))


def along(rng, axis, pl, pt_mean):
    """Momentum with longitudinal component pl along axis and an exponential pT around it."""
    ref = (1.0, 0.0, 0.0) if abs(axis[0]) < 0.9 else (0.0, 1.0, 0.0)
    e1 = unit((axis[1] * ref[2] - axis[2] * ref[1], axis[2] * ref[0] - axis[0] * ref[2],
               axis[0] * ref[1] - axis[1] * ref[0]))
    e2 = (axis[1] * e1[2] - axis[2] * e1[1], axis[2] * e1[0] - axis[0] * e1[2],
          axis[0] * e1[1] - axis[1] * e1[0])
    pt = rng.expovariate(1.0 / pt_mean)
    ph = rng.uniform(0.0, 2 * math.pi)
    return tuple(pl * axis[k] + pt * (math.cos(ph) * e1[k] + math.sin(ph) * e2[k]) for k in range(3))


class Event:
    """Generator-level content of one event."""

    def __init__(self):
        self.mc = []        # dicts: pdg, q, m, p4, vtx, end, status
        self.visible = []   # MC indices of the reconstructable final-state particles
        self.v0s = []       # (pdg, decay vertex, mc index of the two daughters)
        self.pv = (0.0, 0.0, 0.0)
        self.axis = (0.0, 0.0, 1.0)
        self.awkward = None  # kind of awkward event (AWKWARD_KINDS), or None
        self.keep = None     # MC indices of the charged particles to reconstruct (None: all)
        self.special = {}    # MC index -> special treatment of its track (see fill_frame)

    def add(self, name, q, p4, vtx, status=1, sign=1, visible=True, pdg=None, mass=None):
        if pdg is None:
            pdg = PDG[name] * (sign if name not in ("gamma", "nh", "rho", "phi", "K0S") else 1)
        self.mc.append({"pdg": pdg, "q": q, "m": M[name] if mass is None else mass, "p4": p4,
                        "vtx": vtx, "end": (0.0, 0.0, 0.0), "status": status, "name": name})
        idx = len(self.mc) - 1
        if visible and status == 1:
            self.visible.append(idx)
        return idx

    def decay_point(self, rng, idx, ctau):
        part = self.mc[idx]
        p4 = part["p4"]
        p = math.sqrt(p4[0] ** 2 + p4[1] ** 2 + p4[2] ** 2)
        length = rng.expovariate(1.0) * ctau * p / part["m"]
        d = unit(p4)
        end = tuple(part["vtx"][k] + length * d[k] for k in range(3))
        part["end"] = end
        part["status"] = 2
        if idx in self.visible:   # decayed: its daughters are the visible particles
            self.visible.remove(idx)
        return end


def gen_hemisphere(rng, ev, pv, axis, flavour, sign):
    """Particles of one hemisphere; sign = +1 for the quark, -1 for the antiquark side."""
    # leading heavy hadron
    if flavour == 5:
        pB = rng.uniform(0.6, 0.85) * E_BEAM
        b = ev.add("B", 0, p4_of(M["B"], along(rng, axis, pB, 0.3)), pv, sign=sign)
        bv = ev.decay_point(rng, b, CTAU["B"])
        d4, a4 = two_body(rng, ev.mc[b]["p4"], M["Dplus"], M["a1"])
        d = ev.add("Dplus", -sign, d4, bv, sign=-sign)
        a = ev.add("a1", sign, a4, bv, status=2, sign=sign)
        dv = ev.decay_point(rng, d, CTAU["Dplus"])
        ks4, pi4 = two_body(rng, d4, M["Kstar0"], M["pi"])
        ev.add("pi", -sign, pi4, dv, sign=-sign)
        kst = ev.add("Kstar0", 0, ks4, dv, status=2, sign=-sign)
        k4, p4 = two_body(rng, ks4, M["K"], M["pi"])
        ev.add("K", -sign, k4, dv, sign=-sign)
        ev.add("pi", sign, p4, dv, sign=sign)
        ev.mc[kst]["end"] = dv
        r4, pi4 = two_body(rng, a4, M["rho"], M["pi"])
        ev.mc[a]["end"] = bv
        ev.add("pi", sign, pi4, bv, sign=sign)
        r = ev.add("rho", 0, r4, bv, status=2)
        ev.mc[r]["end"] = bv
        x4, y4 = two_body(rng, r4, M["pi"], M["pi"])
        ev.add("pi", 1, x4, bv)
        ev.add("pi", -1, y4, bv, sign=-1)
    elif flavour == 4 or rng.random() < 0.08:
        pc = rng.uniform(0.35, 0.7) * E_BEAM
        if rng.random() < 0.6:   # D*+ -> D0 pi+, D0 -> K- pi+
            ds = ev.add("Dstar", sign, p4_of(M["Dstar"], along(rng, axis, pc, 0.3)), pv, status=2, sign=sign)
            ev.mc[ds]["end"] = pv
            d04, pis4 = two_body(rng, ev.mc[ds]["p4"], M["D0"], M["pi"])
            d0 = ev.add("D0", 0, d04, pv, sign=sign)
            ev.add("pi", sign, pis4, pv, sign=sign)
            dv = ev.decay_point(rng, d0, CTAU["D0"])
            k4, pi4 = two_body(rng, d04, M["K"], M["pi"])
            ev.add("K", -sign, k4, dv, sign=-sign)
            ev.add("pi", sign, pi4, dv, sign=sign)
        else:                    # D+ -> K*0bar pi+, K*0bar -> K- pi+
            d = ev.add("Dplus", sign, p4_of(M["Dplus"], along(rng, axis, pc, 0.3)), pv, sign=sign)
            dv = ev.decay_point(rng, d, CTAU["Dplus"])
            ks4, pi4 = two_body(rng, ev.mc[d]["p4"], M["Kstar0"], M["pi"])
            ev.add("pi", sign, pi4, dv, sign=sign)
            kst = ev.add("Kstar0", 0, ks4, dv, status=2, sign=-sign)
            ev.mc[kst]["end"] = dv
            k4, p4 = two_body(rng, ks4, M["K"], M["pi"])
            ev.add("K", -sign, k4, dv, sign=-sign)
            ev.add("pi", sign, p4, dv, sign=sign)
    # strange hadrons
    if rng.random() < (0.6 if flavour == 3 else 0.35):
        ks = ev.add("K0S", 0, p4_of(M["K0S"], along(rng, axis, rng.expovariate(1 / 4.0) + 0.5, 0.4)), pv)
        kv = ev.decay_point(rng, ks, CTAU["K0S"])
        a4, b4 = two_body(rng, ev.mc[ks]["p4"], M["pi"], M["pi"])
        ia = ev.add("pi", 1, a4, kv)
        ib = ev.add("pi", -1, b4, kv, sign=-1)
        ev.v0s.append((310, kv, ia, ib))
    if rng.random() < 0.15:
        lam = ev.add("Lambda", 0, p4_of(M["Lambda"], along(rng, axis, rng.expovariate(1 / 5.0) + 1.0, 0.4)), pv, sign=sign)
        lv = ev.decay_point(rng, lam, CTAU["Lambda"])
        a4, b4 = two_body(rng, ev.mc[lam]["p4"], M["p"], M["pi"])
        ia = ev.add("p", sign, a4, lv, sign=sign)
        ib = ev.add("pi", -sign, b4, lv, sign=-sign)
        ev.v0s.append((3122, lv, ia, ib))
    if rng.random() < (0.3 if flavour == 3 else 0.12):
        ph = ev.add("phi", 0, p4_of(M["phi"], along(rng, axis, rng.expovariate(1 / 5.0) + 1.0, 0.4)), pv, status=2)
        ev.mc[ph]["end"] = pv
        a4, b4 = two_body(rng, ev.mc[ph]["p4"], M["K"], M["K"])
        ev.add("K", 1, a4, pv)
        ev.add("K", -1, b4, pv, sign=-1)
    # prompt fragmentation: a soft, exponential longitudinal spectrum
    for _ in range(max(2, int(rng.gauss(8.0, 2.5)))):
        u = rng.random()
        name = "pi" if u < 0.8 else ("K" if u < 0.92 else "p")
        q = rng.choice((1, -1))
        p3 = along(rng, axis, rng.expovariate(1 / 1.8) + 0.12, 0.4)
        ev.add(name, q, p4_of(M[name], p3), pv, sign=q)
    if rng.random() < 0.12:
        name = rng.choice(("e", "mu"))
        q = rng.choice((1, -1))
        ev.add(name, q, p4_of(M[name], along(rng, axis, rng.expovariate(1 / 6.0) + 1.0, 0.8)), pv, sign=-q)
    for _ in range(max(1, int(rng.gauss(6.0, 2.0)))):
        ev.add("gamma", 0, p4_of(0.0, along(rng, axis, rng.expovariate(1 / 1.8) + 0.1, 0.25)), pv)
    for _ in range(int(rng.random() < 0.7) + int(rng.random() < 0.2)):
        ev.add("nh", 0, p4_of(M["nh"], along(rng, axis, rng.expovariate(1 / 4.0) + 0.5, 0.4)), pv)


def gen_event(rng, beamspot, hadronic, flavour):
    ev = Event()
    pv = (beamspot[0] + rng.gauss(0.0, 0.010), beamspot[1] + rng.gauss(0.0, 0.0015),
          beamspot[2] + rng.gauss(0.0, 0.8))
    ev.pv = pv
    c = rng.uniform(-0.9, 0.9)
    ph = rng.uniform(0.0, 2 * math.pi)
    axis = (math.sqrt(1 - c * c) * math.cos(ph), math.sqrt(1 - c * c) * math.sin(ph), c)
    ev.axis = axis
    z = ev.add("Z", 0, (0.0, 0.0, 0.0, 2 * E_BEAM), pv, status=110000, visible=False, pdg=23, mass=91.1876)
    ev.mc[z]["end"] = pv
    if not hadronic:            # mu+ mu- pair
        for s in (1, -1):
            ev.add("mu", s, p4_of(M["mu"], tuple(s * E_BEAM * a for a in axis)), pv, sign=-s)
        return ev
    for s in (1, -1):           # the primary quark pair, as the hard-process entries
        i = ev.add("q", s * (2 / 3 if flavour in (2, 4) else -1 / 3),
                   (s * E_BEAM * axis[0], s * E_BEAM * axis[1], s * E_BEAM * axis[2], E_BEAM), pv,
                   status=140001, visible=False, pdg=s * flavour,
                   mass=(0.0, 0.005, 0.002, 0.095, 1.27, 4.18)[flavour])
        ev.mc[i]["end"] = pv
    for s in (1, -1):
        gen_hemisphere(rng, ev, pv, tuple(s * a for a in axis), flavour, s)
    return ev


def make_awkward(rng, ev, kind):
    """Turns a generated hadronic event into an awkward one; the MC truth stays complete."""
    ev.awkward = kind
    if kind in ("one_track", "no_track", "split_singular"):
        for s in (1, -1):         # two hard photons, so that the event still has two jets
            d = unit(tuple(s * a + rng.gauss(0.0, 0.1) for a in ev.axis))
            ev.add("gamma", 0, p4_of(0.0, tuple(rng.uniform(4.0, 12.0) * x for x in d)), ev.pv)
    # prompt charged particles in the barrel, hardest first
    prompt = [i for i in ev.visible if ev.mc[i]["q"] != 0 and ev.mc[i]["vtx"] == ev.pv
              and abs(unit(ev.mc[i]["p4"])[2]) < 0.8]
    prompt.sort(key=lambda i: -(ev.mc[i]["p4"][0] ** 2 + ev.mc[i]["p4"][1] ** 2 + ev.mc[i]["p4"][2] ** 2))
    if kind == "one_track":       # a single track for the vertex fit
        ev.keep = set(prompt[:1])
        ev.special = {i: "clean" for i in prompt[:1]}
    elif kind == "no_track":      # no track at all, or two that fail the track selection
        if rng.random() < 0.5:
            ev.keep = set()
        else:
            ev.keep = set(prompt[:2])
            ev.special = dict(zip(prompt[:2], ("few_tpc_hits", "far_z0")))
    elif kind == "split_singular":
        # P: a track with tiny errors; R: almost collinear with P; Q: roughly perpendicular,
        # from a point 1 mm away in z. The fit of the three converges (Q constrains the
        # vertex along P); Q is pruned, and the refit of P and R is singular along P, where
        # only the beam spot constrains it (below the fit's eigenvalue floor). P runs along
        # x, where the beam spot is widest, through a point 2 mm away from the primary vertex
        # in y, so that its distance to the origin is large (see fill_frame).
        vtx = (ev.pv[0], ev.pv[1] + 0.2 * rng.choice((1, -1)), ev.pv[2])
        phi = rng.uniform(-0.1, 0.1) + (0.0 if rng.random() < 0.5 else math.pi)
        tl = rng.uniform(-0.8, 0.8)
        sp, sr, sq = (rng.choice((1, -1)) for _ in range(3))
        pdir = unit((math.cos(phi), math.sin(phi), tl))
        phr = phi + sr * rng.uniform(0.005, 0.02)
        rdir = unit((math.cos(phr), math.sin(phr), tl + rng.uniform(-0.01, 0.01)))
        phq = phi + sq * (math.pi / 2 + rng.uniform(-0.2, 0.2))
        qdir = unit((math.cos(phq), math.sin(phq), rng.uniform(-0.3, 0.3)))
        qvtx = (vtx[0], vtx[1], vtx[2] + 0.1 * rng.choice((1, -1)))
        ip = ev.add("pi", sp, p4_of(M["pi"], tuple(5.0 * x for x in pdir)), vtx, sign=sp)
        ir = ev.add("pi", -sp, p4_of(M["pi"], tuple(3.0 * x for x in rdir)), vtx, sign=-sp)
        iq = ev.add("pi", sq, p4_of(M["pi"], tuple(8.0 * x for x in qdir)), qvtx, sign=sq)
        ev.keep = {ip, ir, iq}
        ev.special = {ip: "split_precise", ir: "clean", iq: "clean"}
    elif kind == "precise":       # the prompt track farthest from the origin (see fill_frame)
        far = sorted(prompt, key=lambda i: -abs(helix(ev.mc[i]["vtx"], ev.mc[i]["p4"][:3], ev.mc[i]["q"])[0][0]))
        ev.special = {i: kind for i in far[:1]}
    else:                         # duplicate, huge_error: the hardest prompt track
        ev.special = {i: kind for i in prompt[:1]}


# ----------------------------------------------------------------------------------------
# detector response
# ----------------------------------------------------------------------------------------
def helix(vtx, p3, q):
    """Perigee parameters (stored convention) of a particle produced at vtx with momentum p3."""
    pt = math.hypot(p3[0], p3[1])
    kappa = -q * KB / pt
    phiv = math.atan2(p3[1], p3[0])
    cx = vtx[0] - math.sin(phiv) / kappa
    cy = vtx[1] + math.cos(phiv) / kappa
    sgn = 1.0 if kappa > 0 else -1.0
    rc = math.hypot(cx, cy)
    phi0 = math.atan2(-sgn * cx, sgn * cy)
    x0 = cx + math.sin(phi0) / kappa
    y0 = cy - math.cos(phi0) / kappa
    d0 = x0 * math.sin(phi0) - y0 * math.cos(phi0)
    dphi = (phiv - phi0 + math.pi) % (2 * math.pi) - math.pi
    tanl = p3[2] / pt
    z0 = vtx[2] - tanl * dphi / kappa
    return [d0, phi0 % (2 * math.pi), kappa, z0, tanl], rc, 1.0 / abs(kappa)


def track_quality(rng):
    """Hit efficiencies (VDET, ITC, TPC) and the radius where the track leaves the TPC."""
    u = rng.random()
    if u < P_SHORT:
        return (0.9, rng.uniform(0.5, 0.95), rng.uniform(0.6, 0.95), rng.uniform(TPC_RADII[0], 75.0))
    if u < P_SHORT + P_DEGRADED:
        return (rng.uniform(0.6, 0.95), rng.uniform(0.0, 0.9), rng.uniform(0.35, 0.9), 1e9)
    return (0.95, 0.97, 0.95, 1e9)


def hit_counts(rng, vtx, p3, rc, radius, quality):
    eff_v, eff_i, eff_t, r_end = quality
    rv = math.hypot(vtx[0], vtx[1])
    p = math.sqrt(p3[0] ** 2 + p3[1] ** 2 + p3[2] ** 2)
    cos_t = abs(p3[2]) / p
    tanl = p3[2] / math.hypot(p3[0], p3[1])
    rmax = rc + radius
    nvdet = sum(1 for r, acc in VDET_LAYERS if rv < r < rmax and cos_t < acc and rng.random() < eff_v)
    nitc = sum(1 for r in ITC_RADII if rv < r < rmax and cos_t < 0.97 and rng.random() < eff_i)
    ntpc = sum(1 for r in TPC_RADII
               if rv < r < min(rmax, r_end) and abs(vtx[2] + tanl * (r - rv)) < TPC_HALF_Z
               and rng.random() < eff_t)
    return nvdet, nitc, ntpc


def track_covariance(p, pt, kappa, tanl, nvdet, nitc, ntpc):
    """Helix-fit resolutions [cm, rad, 1/cm] of the published size: ~25 um + 95 um/p with the vertex detector."""
    scale = math.sqrt(21.0 / max(ntpc, 4)) * (1.3 if nitc == 0 else 1.0)
    if nvdet >= 1:
        sd0 = math.hypot(0.0025, 0.0095 / p)
        sz0 = math.hypot(0.0028, 0.0100 / p)
    else:
        sd0 = math.hypot(0.015, 0.030 / p)
        sz0 = math.hypot(0.040, 0.060 / p)
    sphi = math.hypot(2e-4, 1e-3 / p) * scale
    stl = math.hypot(5e-4, 1.5e-3 / p) * (1 + tanl * tanl) * scale
    som = abs(kappa) * math.hypot(6e-4 * pt, 0.005) * scale
    sig = np.array([sd0 * scale, sphi, som, sz0 * scale, stl])
    return CORR * np.outer(sig, sig)


def smear(rng, par, cov):
    """Measured parameters: par + a Gaussian residual with covariance cov."""
    return np.array(par) + np.linalg.cholesky(cov) @ np.array([rng.gauss(0, 1) for _ in range(5)])


def dedx_mean(bg):
    b2 = bg * bg / (1.0 + bg * bg)
    return (0.85 + 0.075 * math.log(1.0 + bg * bg)) / b2


# ----------------------------------------------------------------------------------------
# podio output
# ----------------------------------------------------------------------------------------
def _rp_type(rp):
    """ReconstructedParticle type 0, as written by the converter; EDM4hep >= 1.0 has no type."""
    if hasattr(rp, "setType"):
        rp.setType(0)


def fill_frame(rng, ev, mode, run, ievt, classbits, stats, podio, edm4hep, ROOT):
    """All collections of one event, in a podio Frame."""
    tracks = edm4hep.TrackCollection()
    pads = edm4hep.RecDqdxCollection()
    wires = edm4hep.RecDqdxCollection()
    rps = edm4hep.ReconstructedParticleCollection()
    pids = edm4hep.ParticleIDCollection()
    ecal = edm4hep.ClusterCollection()
    hcal = edm4hep.ClusterCollection()
    vertices = edm4hep.VertexCollection()
    jets = edm4hep.ReconstructedParticleCollection()
    v0v = edm4hep.VertexCollection()
    v0p = edm4hep.ReconstructedParticleCollection()
    v0vr = edm4hep.VertexCollection()
    v0pr = edm4hep.ReconstructedParticleCollection()
    refit = edm4hep.ReconstructedParticleCollection()
    is_mc = mode == "mc"

    # MC truth first, so that the links can point at it
    if is_mc:
        mcps = edm4hep.MCParticleCollection()
        mc_objs = []
        for part in ev.mc:
            m = mcps.create()
            m.setPDG(int(part["pdg"]))
            m.setGeneratorStatus(int(part["status"]))
            m.setCharge(float(part["q"]))
            m.setMass(float(part["m"]))
            m.setMomentum(edm4hep.Vector3d(*part["p4"][:3]))
            m.setVertex(edm4hep.Vector3d(*part["vtx"]))
            m.setEndpoint(edm4hep.Vector3d(*part["end"]))
            mc_objs.append(m)
        trk_links = edm4hep.TrackMCParticleLinkCollection()
        cl_links = edm4hep.ClusterMCParticleLinkCollection()

    def add_track(idx, meas, cov, hits, chi2, ndf, dedx_fail):
        """One Track with its state, both dE/dx legs and (MC) its truth link."""
        part = ev.mc[idx]
        t = tracks.create()
        t.setType(0)
        t.setChi2(float(chi2))
        t.setNdf(int(ndf))
        for n in hits:
            t.addToSubdetectorHitNumbers(int(n))
        ts = edm4hep.TrackState()
        ts.location = 0
        ts.D0, ts.phi, ts.omega, ts.Z0, ts.tanLambda = (float(x) for x in meas)
        ts.time = 0.0
        for i in range(5):
            for j in range(i + 1):
                ts.covMatrix[i * (i + 1) // 2 + j] = float(cov[i, j])
        ts.covMatrix[20] = 1e10   # time: not measured
        t.addToTrackStates(ts)
        # dE/dx on both legs; a failed leg stores the track omega (the converter's convention)
        p = math.sqrt(sum(x * x for x in part["p4"][:3]))
        bg = p / part["m"] if part["m"] > 0 else 1e4
        for coll, res, fail in ((pads, 0.10, dedx_fail[0]), (wires, 0.06, dedx_fail[1])):
            d = coll.create()
            q = edm4hep.Quantity()
            sigma = res * math.sqrt(21.0 / max(hits[2], 4))
            if rng.random() < fail:
                q.value = abs(float(meas[2]))
            else:
                q.value = float(dedx_mean(bg) * (1.0 + sigma * rng.gauss(0, 1)))
            q.error = float(dedx_mean(bg) * sigma)
            q.type = 0
            d.setDQdx(q)
            d.setTrack(t)
        if is_mc:
            lk = trk_links.create()
            lk.setFrom(t)
            lk.setTo(mc_objs[idx])
            lk.setWeight(1.0)
        stats["tracks"] += 1
        return t

    # charged particles: tracks, dE/dx, energy-flow objects
    charged = []    # (track object, mc index, measured momentum, ef type)
    for idx in ev.visible:
        part = ev.mc[idx]
        if part["q"] == 0 or (ev.keep is not None and idx not in ev.keep):
            continue
        kind = ev.special.get(idx)    # special treatment in an awkward event, or None
        p3 = part["p4"][:3]
        pt = math.hypot(p3[0], p3[1])
        if pt < 0.05:
            continue
        par, rc, radius = helix(part["vtx"], p3, part["q"])
        nvdet, nitc, ntpc = hit_counts(rng, part["vtx"], p3, rc, radius,
                                       CLEAN if kind else track_quality(rng))
        if kind == "few_tpc_hits":    # fails the TPC-hit cut
            nitc, ntpc = max(nitc, 6), rng.randint(1, 3)
        if ntpc + nitc < 4:
            continue                  # not reconstructed
        p = math.sqrt(sum(x * x for x in p3))
        cov = track_covariance(p, pt, par[2], par[4], nvdet, nitc, ntpc)
        z = [rng.gauss(0, 1) for _ in range(5)]
        tail = 1.0
        if kind is None and rng.random() < P_TAIL:   # mismeasured: residuals the covariance does not describe
            tail = 1.0 + rng.expovariate(1 / 3.0)
            z[0] *= tail
            z[3] *= tail
        meas = np.array(par) + np.linalg.cholesky(cov) @ np.array(z)
        if kind is None:
            u = rng.random()
            if u < P_BAD_COV:         # not positive definite: fails the track selection
                cov[2, 2] = -cov[2, 2]
            elif u < P_BAD_COV + P_FAR_Z0:   # around the |z0| <= 50 cm cut
                meas[3] = rng.choice((-1, 1)) * rng.uniform(40.0, 80.0)
        elif kind == "far_z0":
            meas[3] = rng.choice((-1, 1)) * rng.uniform(55.0, 80.0)
        elif kind in ("precise", "split_precise"):
            # exact parameters with tiny errors: a weight of ~1e21 (~2e16) /cm^2 across the
            # track at the vertex. The phi error is chosen such that phi and d0 move the
            # perigee point by the same amount (along and across the track): the vertex fit
            # projects out the direction along the track, and a much larger weight there
            # would leave a round-off larger than the weight of the other tracks.
            sd = 3e-11 if kind == "precise" else 7e-9
            sa = sd / max(abs(par[0]), 0.002)
            cov = np.diag([sd * sd, sa * sa, sd * sd, sd * sd, sa * sa])
            meas = np.array(par)
        elif kind == "huge_error":
            cov = cov * 1e6
        ndf = max(2 * ntpc + nitc + 2 * nvdet - 5, 1)
        chi2_scale = 1.0 + 0.3 * (tail - 1.0)
        if kind is None and rng.random() < P_CHI2_TAIL:
            chi2_scale *= math.exp(rng.uniform(0.0, math.log(20.0)))
        chi2 = rng.gammavariate(ndf / 2.0, 2.0) * chi2_scale
        t = add_track(idx, meas, cov, (nvdet, nitc, ntpc), chi2, ndf, P_DEDX_FAIL)
        ptm = KB / abs(meas[2])
        pm = tuple(x * ptm / pt for x in p3)
        eftype = {"e": 1, "mu": 2}.get(part["name"], 0)
        charged.append((t, idx, pm, eftype))

        # tracks without an energy-flow object, right after the track of their particle
        if kind == "duplicate":       # nearly identical copy: a collinear pair
            meas2 = meas + 1e-3 * np.sqrt(np.abs(np.diag(cov))) * np.array([rng.gauss(0, 1) for _ in range(5)])
            add_track(idx, meas2, cov, (nvdet, nitc, ntpc), chi2 * rng.uniform(0.9, 1.1), ndf, P_DEDX_FAIL)
            stats["segments"] += 1
        elif kind is None and rng.random() < P_SEGMENT:
            # a short segment of the same particle (a split in the pattern recognition):
            # few TPC hits, poor resolution; mostly fails the track selection
            snv = int(rng.random() < 0.2)
            sni = rng.randint(0, 8)
            snt = min(int(rng.expovariate(1 / 2.0)), 12)
            if snt + sni >= 4:
                scov = track_covariance(p, pt, par[2], par[4], snv, sni, snt) * rng.uniform(2.0, 10.0)
                sndf = max(2 * snt + sni + 2 * snv - 5, 1)
                schi2 = rng.gammavariate(sndf / 2.0, 2.0) * (math.exp(rng.uniform(0.0, math.log(20.0)))
                                                             if rng.random() < 0.3 else 1.5)
                add_track(idx, smear(rng, par, scov), scov, (snv, sni, snt), schi2, sndf, P_DEDX_FAIL_SEGMENT)
                stats["segments"] += 1

    # energy-flow objects: charged ones first, in an order that is NOT the track order,
    # a few tracks without one (as in the converted files), then the neutrals
    order = list(range(len(charged)))
    rng.shuffle(order)
    rp_of_mc = {}
    ef_types = []                 # energy-flow type per RecoParticle, in creation order
    for k in order:
        t, idx, pm, eftype = charged[k]
        if rng.random() < P_NO_EF and idx not in ev.special:
            continue
        mass = {1: M["e"], 2: M["mu"]}.get(eftype, M["pi"])
        r = rps.create()
        _rp_type(r)
        r.setCharge(float(ev.mc[idx]["q"]))
        r.setMass(float(mass))
        r.setMomentum(edm4hep.Vector3f(*pm))
        r.setEnergy(float(math.sqrt(mass ** 2 + sum(x * x for x in pm))))
        r.addToTracks(t)
        rp_of_mc[idx] = (r, eftype)
        ef_types.append(eftype)
    for idx in ev.visible:
        part = ev.mc[idx]
        if part["q"] != 0:
            continue
        is_gamma = part["name"] == "gamma"
        e_true = part["p4"][3]
        e = max(e_true * (1.0 + (0.18 if is_gamma else 0.85) / math.sqrt(e_true) * rng.gauss(0, 1)), 0.05)
        d = unit(part["p4"])
        cl = (ecal if is_gamma else hcal).create()
        cl.setType(0)
        cl.setEnergy(float(e))
        rad = 190.0 if is_gamma else 300.0
        cl.setPosition(edm4hep.Vector3f(*(rad * x for x in d)))
        cl.setITheta(float(math.acos(d[2])))
        cl.setPhi(float(math.atan2(d[1], d[0])))
        r = rps.create()
        _rp_type(r)
        r.setCharge(0.0)
        r.setMass(0.0)
        r.setMomentum(edm4hep.Vector3f(*(e * x for x in d)))
        r.setEnergy(float(e))
        r.addToClusters(cl)
        rp_of_mc[idx] = (r, 4 if is_gamma else 5)
        ef_types.append(4 if is_gamma else 5)
        if is_mc and is_gamma:
            lk = cl_links.create()
            lk.setFrom(cl)
            lk.setTo(mc_objs[idx])
            lk.setWeight(1.0)
    # ParticleID index-parallel to RecoParticles
    for r, eftype in zip(rps, ef_types):
        pid = pids.create()
        pid.setType(int(eftype))
        pid.setPDG(0)
        pid.setAlgorithmType(0)
        pid.setLikelihood(0.0)
        pid.setParticle(r)

    # stored primary vertex and a two-jet split along the event axis
    v = vertices.create()
    v.setType(2)
    v.setPosition(edm4hep.Vector3f(*(ev.pv[k] + rng.gauss(0.0, 0.002) for k in range(3))))
    v.setCovMatrix(edm4hep.CovMatrix3f(4e-6, 0.0, 4e-6, 0.0, 0.0, 4e-6))
    for idx, (r, _) in rp_of_mc.items():
        if ev.mc[idx]["vtx"] == ev.pv and ev.mc[idx]["q"] != 0:
            v.addToParticles(r)
    axis = None
    for part in ev.mc:
        if part["status"] == 140001:
            axis = unit(part["p4"])
            break
    if axis is not None:
        for s in (1, -1):
            j = jets.create()
            mom = [0.0, 0.0, 0.0]
            en = 0.0
            for idx, (r, _) in rp_of_mc.items():
                pm = ev.mc[idx]["p4"]
                if s * (pm[0] * axis[0] + pm[1] * axis[1] + pm[2] * axis[2]) > 0:
                    m3 = r.getMomentum()
                    mom = [mom[0] + m3.x, mom[1] + m3.y, mom[2] + m3.z]
                    en += r.getEnergy()
                    j.addToParticles(r)
            j.setMomentum(edm4hep.Vector3f(*mom))
            j.setEnergy(float(en))

    # V0 candidates of the converter: the generated K0s/Lambda with both daughters measured
    trk_of_mc = {idx: t for t, idx, _, _ in charged}
    for pdg, dv, ia, ib in ev.v0s:
        if ia not in trk_of_mc or ib not in trk_of_mc:
            continue
        for vcoll, pcoll in ((v0v, v0p), (v0vr, v0pr)):
            vv = vcoll.create()
            vv.setType(0)
            vv.setPosition(edm4hep.Vector3f(*(dv[k] + rng.gauss(0.0, 0.05) for k in range(3))))
            vp = pcoll.create()
            vp.setPDG(pdg)
            _rp_type(vp)
            p3 = [ev.mc[ia]["p4"][k] + ev.mc[ib]["p4"][k] for k in range(3)]
            vp.setMomentum(edm4hep.Vector3f(*p3))
            vp.setMass(float(M["K0S"] if pdg == 310 else M["Lambda"]))
            vp.setEnergy(float(ev.mc[ia]["p4"][3] + ev.mc[ib]["p4"][3]))
            vp.addToTracks(trk_of_mc[ia])
            vp.addToTracks(trk_of_mc[ib])
            vp.setDecayVertex(vv)
            vv.addToParticles(vp)
        for i in (ia, ib):
            rf = refit.create()
            _rp_type(rf)
            rf.setCharge(float(ev.mc[i]["q"]))
            rf.setMomentum(edm4hep.Vector3f(*ev.mc[i]["p4"][:3]))
            rf.setEnergy(float(ev.mc[i]["p4"][3]))
            rf.addToTracks(trk_of_mc[i])

    header = edm4hep.EventHeaderCollection()
    h = header.create()
    h.setEventNumber(int(ievt))
    h.setRunNumber(int(run))
    h.setTimeStamp(int(1000 + 10 * ievt))
    h.setWeight(0.0)              # not filled by the converter either
    cb = ROOT.podio.UserDataCollection["uint32_t"]()
    cb.push_back(int(classbits))
    runinfo = ROOT.podio.UserDataCollection["float"]()
    for x in (E_BEAM * 1e3, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0):   # placeholders, not read by stage1
        runinfo.push_back(float(x))

    frame = podio.Frame()
    for coll, name in ((cb, "ClassBitset"), (pads, "dEdxPads"), (wires, "dEdxWires"),
                       (ecal, "ElectromagneticClusters"), (header, "EventHeader"),
                       (hcal, "HadronicClusters"), (jets, "Jets"), (pids, "ParticleID"),
                       (rps, "RecoParticles"), (refit, "RefittedParticles"),
                       (runinfo, "RunInformation"), (tracks, "Tracks"), (v0p, "V0Particles"),
                       (v0pr, "V0ParticlesRefitted"), (v0v, "V0Vertices"),
                       (v0vr, "V0VerticesRefitted"), (vertices, "Vertices")):
        frame.put(coll, name)
    if is_mc:
        gen = edm4hep.GeneratorEventParametersCollection()
        g = gen.create()
        g.setSqrts(2 * E_BEAM)
        frame.put(cl_links, "ecalClusterMCLink")
        frame.put(gen, "GeneratorInfo")
        frame.put(mcps, "MCParticles")
        frame.put(trk_links, "trackMCLink")
    return frame


def class_bits(rng, hadronic):
    bits = 0
    for b in (0, 1, 2, 3, 5, 10, 11, 12, 13, 16, 17):
        if rng.random() < 0.5:
            bits |= 1 << b
    return bits | (1 << (CLASS_HADRONIC_BIT if hadronic else CLASS_LEPTONPAIR_BIT))


def new_stats(runs):
    s = dict(events=0, hadronic=0, bb=0, tracks=0, segments=0)
    s.update({f"run {r}": 0 for r in runs})
    s.update({f"awkward {k}": 0 for k in AWKWARD_KINDS})
    return s


def write_chunk(task):
    """Events [start, stop) of one file into the podio file `path`; returns the counts."""
    path, mode, start, stop, seed, runs, beamspots = task
    import ROOT
    import podio
    import edm4hep
    from podio.root_io import Writer

    stats = new_stats(runs)
    writer = Writer(path)
    for i in range(start, stop):
        rng = random.Random(f"{seed}/{mode}/{i}")      # per-event stream: event i depends on nothing else
        hadronic = rng.random() > P_MUMU
        u = rng.random()
        flavour = 1 + sum(u >= c for c in FLAVOUR_CUMULATIVE)   # 1 d, 2 u, 3 s, 4 c, 5 b
        u, k = rng.random(), 0
        while k < len(runs) - 1 and u >= sum(DATA_RUN_FRACTIONS[:k + 1]):
            k += 1
        ev = gen_event(rng, beamspots[k], hadronic, flavour)
        if hadronic and rng.random() < AWKWARD_FRACTION:
            make_awkward(rng, ev, AWKWARD_KINDS[int(rng.random() * len(AWKWARD_KINDS))])
            stats[f"awkward {ev.awkward}"] += 1
        frame = fill_frame(rng, ev, mode, runs[k], i + 1, class_bits(rng, hadronic), stats, podio, edm4hep, ROOT)
        stats["events"] += 1
        stats["hadronic"] += int(hadronic)
        stats["bb"] += int(hadronic and flavour == 5)
        stats[f"run {runs[k]}"] += 1
        writer.write_frame(frame, "events")
    writer._writer.finish()
    del writer._writer    # podio's atexit hook would finish it a second time
    return stats


def join_chunks(task):
    """Joins the chunk files, in order, into one podio file (events tree + metadata of the first)."""
    paths, out = task
    import ROOT
    chain = ROOT.TChain("events")
    for p in paths:
        chain.Add(p)
    fout = ROOT.TFile(out, "RECREATE")
    if chain.Merge(fout, 0, "fast keep") <= 0:
        raise RuntimeError(f"could not join {len(paths)} chunk files into {out}")
    fin = ROOT.TFile.Open(paths[0])
    meta = fin.Get("podio_metadata")
    fout.cd()
    meta.CloneTree(-1, "fast").Write()
    fout.Close()
    fin.Close()
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", required=True, help="output directory")
    ap.add_argument("--nevents", type=int, required=True, help="events per file")
    ap.add_argument("--seed", default="1", help="random seed (any string)")
    ap.add_argument("--jobs", type=int, default=len(os.sched_getaffinity(0)),
                    help="worker processes (the files do not depend on it; default: available CPUs)")
    args = ap.parse_args()
    t0 = time.time()

    # data runs: three selected runs spread over the year and one rejected run of the run list
    sys.path.insert(0, os.path.join(REPO, "src"))
    import run_list
    good = sorted(run_list.good_runs(DATA_YEAR))
    bad = sorted(set(run_list.load(DATA_YEAR)) - set(good))
    if len(good) < 3 or not bad:
        sys.exit("make_synthetic_input: the run list needs three selected runs and a rejected run")
    data_runs = (good[0], good[len(good) // 2], good[-1], bad[0])
    os.makedirs(args.out, exist_ok=True)
    beamspot_file = os.path.join(args.out, "beamspot.json")
    with open(beamspot_file, "w") as fh:
        json.dump({str(r): dict(zip("xyz", bs)) for r, bs in zip(data_runs, DATA_BEAMSPOTS)}, fh, indent=1)

    # chunks of CHUNK events, written by the workers and joined in order
    tmp = os.path.join(args.out, "chunks.tmp")
    shutil.rmtree(tmp, ignore_errors=True)
    os.makedirs(tmp)
    files = {"data": os.path.join(args.out, DATA_YEAR, "synthetic_data.root"),
             "mc": os.path.join(args.out, "QQB", "synthetic_mc.root")}
    setup = {"data": (data_runs, DATA_BEAMSPOTS), "mc": ((1,), ((0.0, 0.0, 0.0),))}
    tasks, chunks = [], {}
    for mode in files:
        os.makedirs(os.path.dirname(files[mode]), exist_ok=True)
        for start in range(0, args.nevents, CHUNK):
            path = os.path.join(tmp, f"{mode}_{start // CHUNK:05d}.root")
            tasks.append((path, mode, start, min(start + CHUNK, args.nevents), args.seed, *setup[mode]))
            chunks.setdefault(mode, []).append(path)
    # PyROOT compiles its bindings at first use (~10 s): once here, for both modes, so that
    # the forked workers inherit them (a throwaway event, with its own random stream)
    import ROOT
    import podio
    import edm4hep
    for mode in files:
        rng = random.Random("warm-up")
        fill_frame(rng, gen_event(rng, (0.0, 0.0, 0.0), True, 5), mode, 1, 1, 0, new_stats(()),
                   podio, edm4hep, ROOT)
    jobs = max(1, min(args.jobs, len(tasks)))
    if jobs == 1:
        results = [write_chunk(t) for t in tasks]
        for mode in files:
            join_chunks((chunks[mode], files[mode]))
    else:
        with multiprocessing.get_context("fork").Pool(jobs) as pool:
            results = pool.map(write_chunk, tasks, chunksize=1)
            pool.map(join_chunks, [(chunks[mode], files[mode]) for mode in files])
    shutil.rmtree(tmp)

    stats = {mode: new_stats(setup[mode][0]) for mode in files}
    for task, res in zip(tasks, results):
        for key, n in res.items():
            stats[task[1]][key] += n
    for mode in files:
        s = stats[mode]
        runs = ", ".join(f"{s[f'run {r}']} in run {r}" for r in setup[mode][0]) if mode == "data" else ""
        awk = sum(s[f"awkward {k}"] for k in AWKWARD_KINDS)
        print(f"synthetic {mode + ':':5s} {s['events']} events, {s['hadronic']} hadronic, {s['bb']} bbbar, "
              f"{awk} awkward, {s['tracks'] / max(s['events'], 1):.1f} tracks/event "
              f"({s['segments'] / max(s['events'], 1):.1f} without an energy-flow object)")
        if runs:
            print(f"  runs: {runs} (the last one rejected by the run list)")
        print("  awkward: " + ", ".join(f"{k} {s[f'awkward {k}']}" for k in AWKWARD_KINDS))
    print(f"beam-spot file: {beamspot_file}")
    print(f"written in {time.time() - t0:.1f} s with {jobs} worker process(es)")


if __name__ == "__main__":
    main()
