#!/usr/bin/env python3
"""
compare_hadrons.py

Heavy-flavour hadron content comparison of three samples:
  Paper   : ALEPH 1994 qqbar MC           (event flavour from the MC record)
  HVFL    : HVFL 2026 production          (event flavour from the MC record)
  Pythia8 : py8Zff Zbb/Zcc/Zss/Zuu/Zdd    (event flavour = sample; record flavour cross-checked)

Replaces extract_hadrons.cpp + md.py and fixes:
  * every species is counted together with its charge conjugate
  * generator copies and B0/Bs0 oscillation entries are counted once
    (a particle with a daughter of the same |PDG| is skipped)
  * event flavour = quark whose parent is the Z (fallback: first quark in record)
  * exact <n> per event (no histogram overflow clipping)
  * errors from the event-by-event spread, sigma = sqrt(Var(n)/N),
    computed directly for species, categories and totals (correlations included)
  * ratios to Paper with errors; |pull| > 3 highlighted
  * sanity check: weakly-decaying b (c) hadrons per b (c) event should be ~2

Outputs in --outdir:
  comparison.md        all tables
  multiplicities.csv   every observable (incl. particle / antiparticle separately)
  plots/*.png          per flavour and category: <n> + ratio-to-Paper panel
  raw_counts.npz       cache of the event loop (reused unless --reprocess)

Run in the same key4hep / LCG environment used to build extract_hadrons:
  python3 compare_hadrons.py --outdir tables --threads 8
  python3 compare_hadrons.py --outdir test --max-files 3      # quick test
"""

import argparse
import csv
import glob
import json
import os
import sys
from pathlib import Path

import numpy as np

# ----------------------------------------------------------------------------
#  Configuration
# ----------------------------------------------------------------------------
PAPER_DIR = "/eos/experiment/aleph/EDM4HEP/MC/1994/QQB"
HVFL_DIR = "/eos/experiment/fcc/ee/analyses/case-studies/aleph/FULLSIM/QQB"    #"/eos/user/h/hfatehi/2026"
PY8_BASE = "/eos/experiment/fcc/ee/analyses/case-studies/aleph/FULLSIM/py8Zff/edm4hep"
PY8_SUBDIR = {5: "Zbb", 4: "Zcc", 3: "Zss", 2: "Zuu", 1: "Zdd"}
FILE_PREFIX = "file://"

SAMPLES = ["Paper", "HVFL", "Pythia8"]
REF = "Paper"
FLAV_ORDER = [5, 4, 3, 2, 1]
FLAV_NAME = {5: "b", 4: "c", 3: "s", 2: "u", 1: "d"}
R_SM = {5: 0.216, 4: 0.172, 3: 0.220, 2: 0.172, 1: 0.220}   # Z -> qq fractions

CATEGORY_ORDER = ["B_meson", "D_meson", "b_baryon", "c_baryon"]
CAT_TITLE = {"B_meson": "B mesons", "D_meson": "D mesons",
             "b_baryon": "b baryons", "c_baryon": "c baryons"}

# particle only; the charge conjugate is added automatically
SPECIES = {
    "B_meson": [("B0", 511), ("B+", 521), ("Bs0", 531), ("Bc+", 541),
                ("B*0", 513), ("B*+", 523), ("Bs*0", 533), ("Bc*+", 543)],
    "D_meson": [("D+", 411), ("D0", 421), ("Ds+", 431),
                ("D*+", 413), ("D*0", 423), ("Ds*+", 433),
                ("D1(2420)+", 10413), ("D1(2420)0", 10423),
                ("D2*(2460)+", 415), ("D2*(2460)0", 425)],
    "b_baryon": [("Lambda_b0", 5122),
                 ("Sigma_b-", 5112), ("Sigma_b0", 5212), ("Sigma_b+", 5222),
                 ("Sigma_b*-", 5114), ("Sigma_b*0", 5214), ("Sigma_b*+", 5224),
                 ("Xi_b-", 5132), ("Xi_b0", 5232),
                 ("Xi_b'-", 5312), ("Xi_b'0", 5322),
                 ("Xi_b*-", 5314), ("Xi_b*0", 5324),
                 ("Omega_b-", 5332), ("Omega_b*-", 5334)],
    "c_baryon": [("Lambda_c+", 4122),
                 ("Sigma_c0", 4112), ("Sigma_c+", 4212), ("Sigma_c++", 4222),
                 ("Sigma_c*0", 4114), ("Sigma_c*+", 4214), ("Sigma_c*++", 4224),
                 ("Xi_c0", 4132), ("Xi_c+", 4232),
                 ("Xi_c'0", 4312), ("Xi_c'+", 4322),
                 ("Xi_c*0", 4314), ("Xi_c*+", 4324),
                 ("Omega_c0", 4332), ("Omega_c*0", 4334)],
}
WEAK_B = [511, 521, 531, 541, 5122, 5132, 5232, 5332]
WEAK_C = [411, 421, 431, 4122, 4132, 4232, 4332]

# ----------------------------------------------------------------------------
#  C++ event loop (JIT-compiled through PyROOT)
# ----------------------------------------------------------------------------
CPP = r"""
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <edm4hep/MCParticleData.h>
#include <podio/ObjectID.h>
#include <cmath>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace hc {
using ROOT::VecOps::RVec;
using MCVec = RVec<edm4hep::MCParticleData>;
using IdVec = RVec<podio::ObjectID>;

struct Result {
  std::vector<double> nEv, s1, s2;       // nEv[6], s1/s2[6*nObs]  (index = flavour code)
  double nTotal = 0, nSel = 0, nNoFlav = 0, nFallback = 0, nAgree = 0, nDisagree = 0;
  std::string info;
};

// +a : quark with a Z parent found;  -a : fallback to first quark;  0 : no quark
inline int recordFlavour(const MCVec& mc, const IdVec& par) {
  int first = 0;
  const int n = mc.size();
  for (int i = 0; i < n; ++i) {
    const int a = std::abs(mc[i].PDG);
    if (a < 1 || a > 5) continue;
    if (first == 0) first = a;
    for (unsigned j = mc[i].parents_begin; j < mc[i].parents_end && j < par.size(); ++j) {
      const int p = par[j].index;
      if (p >= 0 && p < n && mc[p].PDG == 23) return a;
    }
  }
  return -first;
}

// copy or oscillation entry: has a daughter with the same |PDG|
inline bool isCopy(const MCVec& mc, const IdVec& dau, int i) {
  const int a = std::abs(mc[i].PDG);
  for (unsigned j = mc[i].daughters_begin; j < mc[i].daughters_end && j < dau.size(); ++j) {
    const int d = dau[j].index;
    if (d >= 0 && d < (int)mc.size() && std::abs(mc[d].PDG) == a) return true;
  }
  return false;
}

Result run(const std::vector<std::string>& files, int fixedFlavour,
           const std::vector<int>& pairPdg, const std::vector<int>& pairObs, int nObs)
{
  std::unordered_map<int, std::vector<int>> pdgToObs;
  for (size_t k = 0; k < pairPdg.size(); ++k) pdgToObs[pairPdg[k]].push_back(pairObs[k]);

  ROOT::RDataFrame df0("events", files);
  ROOT::RDF::RNode df = df0;
  Result R;

  auto pick = [&](const std::vector<std::string>& cands) {
    for (const auto& c : cands) if (df0.HasColumn(c)) return c;
    return std::string();
  };
  std::string par  = pick({"_MCParticles_parents", "MCParticles#0"});
  std::string dau  = pick({"_MCParticles_daughters", "MCParticles#1"});
  std::string bits = df0.HasColumn("ClassBitset") ? std::string("ClassBitset") : std::string();

  R.info = "parents=" + (par.empty() ? std::string("MISSING") : par) +
           " daughters=" + (dau.empty() ? std::string("MISSING") : dau) +
           " ClassBitset=" + (bits.empty() ? std::string("MISSING(no selection!)") : bits);

  if (par.empty())  { df = df.Define("hc_par",  [] { return IdVec{}; });                     par  = "hc_par"; }
  if (dau.empty())  { df = df.Define("hc_dau",  [] { return IdVec{}; });                     dau  = "hc_dau"; }
  if (bits.empty()) { df = df.Define("hc_bits", [] { return RVec<uint32_t>{0xFFFFFFFFu}; }); bits = "hc_bits"; }

  struct Acc {
    std::vector<double> nEv, s1, s2; std::vector<int> cnt;
    double tot = 0, sel = 0, noflav = 0, fb = 0, agree = 0, dis = 0;
  };
  const unsigned nSlots = df0.GetNSlots();
  std::vector<Acc> acc(nSlots);
  for (auto& A : acc) {
    A.nEv.assign(6, 0.); A.s1.assign(6 * nObs, 0.); A.s2.assign(6 * nObs, 0.); A.cnt.assign(nObs, 0);
  }

  df.ForeachSlot(
    [&](unsigned int slot, const RVec<uint32_t>& b, const MCVec& mc,
        const IdVec& pr, const IdVec& dt)
    {
      Acc& A = acc[slot];
      A.tot += 1;
      if (b.empty() || !(b[0] & (1u << 15))) return;
      A.sel += 1;

      const int code = recordFlavour(mc, pr);
      int fl;
      if (fixedFlavour > 0) {
        fl = fixedFlavour;
        if (std::abs(code) == fixedFlavour) A.agree += 1; else A.dis += 1;
        if (code < 0) A.fb += 1;
      } else {
        if (code == 0) { A.noflav += 1; return; }
        if (code < 0) A.fb += 1;
        fl = std::abs(code);
      }
      A.nEv[fl] += 1;

      std::fill(A.cnt.begin(), A.cnt.end(), 0);
      for (int i = 0; i < (int)mc.size(); ++i) {
        auto it = pdgToObs.find(mc[i].PDG);
        if (it == pdgToObs.end()) continue;
        if (isCopy(mc, dt, i)) continue;
        for (int o : it->second) A.cnt[o] += 1;
      }
      double* s1 = &A.s1[fl * nObs];
      double* s2 = &A.s2[fl * nObs];
      for (int o = 0; o < nObs; ++o) { const double c = A.cnt[o]; s1[o] += c; s2[o] += c * c; }
    },
    {bits, "MCParticles", par, dau});

  R.nEv.assign(6, 0.); R.s1.assign(6 * nObs, 0.); R.s2.assign(6 * nObs, 0.);
  for (const auto& A : acc) {
    for (int k = 0; k < 6; ++k) R.nEv[k] += A.nEv[k];
    for (size_t k = 0; k < R.s1.size(); ++k) { R.s1[k] += A.s1[k]; R.s2[k] += A.s2[k]; }
    R.nTotal += A.tot; R.nSel += A.sel; R.nNoFlav += A.noflav;
    R.nFallback += A.fb; R.nAgree += A.agree; R.nDisagree += A.dis;
  }
  return R;
}
}  // namespace hc
"""


# ----------------------------------------------------------------------------
#  Observables
# ----------------------------------------------------------------------------
def build_observables():
    obs = []

    def add(**kw):
        kw["i"] = len(obs)
        obs.append(kw)

    for cat in CATEGORY_ORDER:
        allp = []
        for name, pdg in SPECIES[cat]:
            add(key=f"{name} + c.c.", kind="cc", cat=cat, species=name, pdg=pdg, pdgs=[pdg, -pdg])
            add(key=name, kind="particle", cat=cat, species=name, pdg=pdg, pdgs=[pdg])
            add(key=f"anti-{name}", kind="antiparticle", cat=cat, species=name, pdg=-pdg, pdgs=[-pdg])
            allp += [pdg, -pdg]
        add(key=f"all {CAT_TITLE[cat]}", kind="category", cat=cat, species="", pdg=0, pdgs=allp)
    add(key="weakly-decaying b hadrons", kind="weak", cat="", species="", pdg=0,
        pdgs=[s * p for p in WEAK_B for s in (1, -1)])
    add(key="weakly-decaying c hadrons", kind="weak", cat="", species="", pdg=0,
        pdgs=[s * p for p in WEAK_C for s in (1, -1)])
    return obs


# ----------------------------------------------------------------------------
#  ROOT part
# ----------------------------------------------------------------------------
def setup_root(threads):
    import ROOT
    ROOT.gROOT.SetBatch(True)
    old = ROOT.gErrorIgnoreLevel
    ROOT.gErrorIgnoreLevel = ROOT.kFatal
    for lib in ("libpodio", "libpodioDict", "libedm4hep", "libedm4hepDict"):
        ROOT.gSystem.Load(lib)
    ROOT.gErrorIgnoreLevel = old
    if threads > 1:
        ROOT.EnableImplicitMT(threads)
    if not ROOT.gInterpreter.Declare(CPP):
        sys.exit("C++ helper failed to compile - is the key4hep/edm4hep environment sourced?")
    return ROOT


def scan_files(ROOT, d, max_files):
    files = sorted(glob.glob(os.path.join(d, "*.root")))
    if max_files:
        files = files[:max_files]
    ok, bad = [], []
    old = ROOT.gErrorIgnoreLevel
    ROOT.gErrorIgnoreLevel = ROOT.kFatal
    for f in files:
        tf = ROOT.TFile.Open(FILE_PREFIX + f)
        good = bool(tf) and not tf.IsZombie()
        if good:
            t = tf.Get("events")
            good = bool(t) and t.GetEntries() > 0
        if tf:
            tf.Close()
        (ok if good else bad).append(f)
    ROOT.gErrorIgnoreLevel = old
    return ok, bad


def run_event_loops(args, obs):
    ROOT = setup_root(args.threads)
    nObs = len(obs)
    pdgs, oidx = ROOT.std.vector["int"](), ROOT.std.vector["int"]()
    for o in obs:
        for p in o["pdgs"]:
            pdgs.push_back(p)
            oidx.push_back(o["i"])

    R = {s: dict(nEv=np.zeros(6), s1=np.zeros((6, nObs)), s2=np.zeros((6, nObs))) for s in SAMPLES}
    runs = []
    jobs = [("Paper", PAPER_DIR, 0), ("HVFL", HVFL_DIR, 0)] + \
           [("Pythia8", os.path.join(PY8_BASE, PY8_SUBDIR[f]), f) for f in FLAV_ORDER]

    for sample, d, fixed in jobs:
        tag = sample if fixed == 0 else f"Pythia8 ({FLAV_NAME[fixed]})"
        info = dict(tag=tag, sample=sample, dir=d, fixed=fixed, files_ok=0, files_bad=0,
                    total=0, sel=0, noflav=0, fallback=0, agree=0, disagree=0,
                    columns="", status="ok")
        runs.append(info)
        if not os.path.isdir(d):
            info["status"] = "missing directory"
            print(f"[{tag}] missing directory {d}", file=sys.stderr)
            continue
        ok, bad = scan_files(ROOT, d, args.max_files)
        info["files_ok"], info["files_bad"] = len(ok), len(bad)
        for b in bad:
            print(f"[{tag}] skipping unreadable/empty file {b}", file=sys.stderr)
        if not ok:
            info["status"] = "no readable files"
            continue

        files = ROOT.std.vector["std::string"]()
        for f in ok:
            files.push_back(FILE_PREFIX + f)
        print(f"[{tag}] processing {len(ok)} files ...", file=sys.stderr)
        try:
            res = ROOT.hc.run(files, fixed, pdgs, oidx, nObs)
        except Exception as e:  # e.g. unexpected column type
            info["status"] = f"FAILED: {e}"
            print(f"[{tag}] FAILED: {e}", file=sys.stderr)
            continue

        R[sample]["nEv"] += np.array(list(res.nEv), dtype=float)
        R[sample]["s1"] += np.array(list(res.s1), dtype=float).reshape(6, nObs)
        R[sample]["s2"] += np.array(list(res.s2), dtype=float).reshape(6, nObs)
        info.update(total=res.nTotal, sel=res.nSel, noflav=res.nNoFlav,
                    fallback=res.nFallback, agree=res.nAgree, disagree=res.nDisagree,
                    columns=str(res.info))
        print(f"[{tag}] {int(res.nSel):,} selected events; {res.info}", file=sys.stderr)
    return R, runs


# ----------------------------------------------------------------------------
#  Cache
# ----------------------------------------------------------------------------
def save_cache(path, R, runs, obs, max_files):
    arrs = {f"{s}__{k}": R[s][k] for s in SAMPLES for k in ("nEv", "s1", "s2")}
    meta = json.dumps(dict(runs=runs, max_files=max_files))
    np.savez(path, obs_keys=np.array([o["key"] for o in obs]), meta=np.array(meta), **arrs)


def load_cache(path, obs, max_files):
    d = np.load(path)
    if list(d["obs_keys"]) != [o["key"] for o in obs]:
        return None
    meta = json.loads(str(d["meta"]))
    if meta.get("max_files") != max_files:
        return None
    R = {s: {k: d[f"{s}__{k}"] for k in ("nEv", "s1", "s2")} for s in SAMPLES}
    return R, meta["runs"]


# ----------------------------------------------------------------------------
#  Statistics
# ----------------------------------------------------------------------------
def compute(R):
    S = {}
    for s in SAMPLES:
        N = R[s]["nEv"][:, None]
        with np.errstate(divide="ignore", invalid="ignore"):
            mean = R[s]["s1"] / N
            var = R[s]["s2"] / N - mean ** 2
            err = np.sqrt(np.clip(var, 0, None) / N)
        mean[R[s]["nEv"] == 0, :] = np.nan
        err[R[s]["nEv"] == 0, :] = np.nan
        S[s] = dict(N=R[s]["nEv"], mean=mean, err=err)
    return S


def fmt(m, e, scale=1.0):
    if not np.isfinite(m):
        return "—"
    m, e = m * scale, e * scale
    if m == 0:
        return "0"
    if not np.isfinite(e) or e <= 0:
        return f"{m:.3g}"
    d = int(min(6, max(0, 1 - np.floor(np.log10(e)))))
    return f"{m:.{d}f} ± {e:.{d}f}"


def ratio(a, ea, b, eb):
    """returns (r, er, pull) or None"""
    if not (np.isfinite(a) and np.isfinite(b)) or b <= 0:
        return None
    r = a / b
    er = r * np.hypot(ea / a, eb / b) if a > 0 else 0.0
    den = np.hypot(ea, eb)
    pull = (a - b) / den if den > 0 else 0.0
    return r, er, pull


def ratio_cell(a, ea, b, eb):
    x = ratio(a, ea, b, eb)
    if x is None:
        return "—"
    r, er, pull = x
    if r == 0:
        return "0"
    txt = f"{r:.3f} ± {er:.3f}"
    return f"**{txt}**" if abs(pull) > 3 else txt


# ----------------------------------------------------------------------------
#  Markdown
# ----------------------------------------------------------------------------
def md_table(header, rows):
    out = ["| " + " | ".join(header) + " |", "|" + "|".join(["---"] * len(header)) + "|"]
    out += ["| " + " | ".join(r) + " |" for r in rows]
    return out


def write_markdown(path, S, obs, runs):
    L = ["# Heavy-flavour hadron content: Paper (LEP1994) vs HVFL2026 vs Pythia8", ""]
    L += ["Values are the mean number of hadrons per selected hadronic event "
          "(ClassBitset bit 15), both hemispheres, particle + antiparticle. "
          "Generator copies and B⁰/B_s⁰ oscillation entries are counted once. "
          "Uncertainties are statistical, σ = √(Var(n)/N). Ratios are to Paper; "
          "**bold** marks |pull| > 3.", ""]

    # ---- samples overview
    L += ["## Samples", ""]
    rows = []
    for r in runs:
        sel = r["sel"] or 0
        fb = f"{r['fallback'] / sel:.2%}" if sel else "—"
        nf = f"{r['noflav'] / sel:.2%}" if sel and r["fixed"] == 0 else "—"
        n_ad = r["agree"] + r["disagree"]
        dis = f"{r['disagree'] / n_ad:.2%}" if n_ad else "—"
        rows.append([r["tag"], str(r["files_ok"]), str(r["files_bad"]),
                     f"{int(r['total']):,}", f"{int(sel):,}", nf, fb, dis,
                     r["status"], f"`{r['columns']}`" if r["columns"] else ""])
    L += md_table(["Sample", "Files OK", "Files bad", "Events read", "Pass bit 15",
                   "No flavour", "Flavour by fallback", "Record ≠ nominal", "Status", "Columns"], rows)
    L += ["", "*Flavour by fallback*: no quark with a Z parent was found, the first quark "
          "in the record was used. *Record ≠ nominal*: Pythia8 events whose record flavour "
          "differs from the sample flavour.", ""]

    # ---- flavour composition
    L += ["## Flavour composition", ""]
    rows = []
    tot = {s: S[s]["N"][1:6].sum() for s in SAMPLES}
    for fl in FLAV_ORDER:
        row = [FLAV_NAME[fl], f"{R_SM[fl]:.3f}"]
        for s in ("Paper", "HVFL"):
            row.append(f"{S[s]['N'][fl] / tot[s]:.4f}" if tot[s] else "—")
        row += [f"{int(S[s]['N'][fl]):,}" for s in SAMPLES]
        rows.append(row)
    L += md_table(["Flavour", "SM R_q", "Paper fraction", "HVFL fraction",
                   "Paper events", "HVFL events", "Pythia8 events"], rows)
    L.append("")

    cat_rows = [o for o in obs if o["kind"] in ("category", "weak")]
    for fl in FLAV_ORDER:
        L += [f"## {FLAV_NAME[fl]}-jets", ""]
        L.append(" | ".join(f"**{s}**: {int(S[s]['N'][fl]):,} evt" for s in SAMPLES))
        L.append("")

        # categories (per event, not scaled)
        L += ["### Category totals (hadrons per event)", ""]
        rows = []
        for o in cat_rows:
            i = o["i"]
            m = {s: S[s]["mean"][fl, i] for s in SAMPLES}
            e = {s: S[s]["err"][fl, i] for s in SAMPLES}
            rows.append([o["key"]] + [fmt(m[s], e[s]) for s in SAMPLES] +
                        [ratio_cell(m[s], e[s], m[REF], e[REF]) for s in SAMPLES if s != REF])
        L += md_table(["Category"] + [f"{s} ⟨n⟩" for s in SAMPLES] +
                      [f"{s} / {REF}" for s in SAMPLES if s != REF], rows)
        L.append("")

        # species
        for cat in CATEGORY_ORDER:
            rows = []
            for o in obs:
                if o["kind"] != "cc" or o["cat"] != cat:
                    continue
                i = o["i"]
                m = {s: S[s]["mean"][fl, i] for s in SAMPLES}
                e = {s: S[s]["err"][fl, i] for s in SAMPLES}
                if all((not np.isfinite(m[s])) or m[s] == 0 for s in SAMPLES):
                    continue
                rows.append([f"`{o['species']}`", str(o["pdg"])] +
                            [fmt(m[s], e[s], 1e3) for s in SAMPLES] +
                            [ratio_cell(m[s], e[s], m[REF], e[REF]) for s in SAMPLES if s != REF])
            if not rows:
                continue
            L += [f"### {CAT_TITLE[cat]} (+ c.c.)", ""]
            L += md_table(["Species", "PDG"] + [f"{s} ⟨n⟩ ×10³" for s in SAMPLES] +
                          [f"{s} / {REF}" for s in SAMPLES if s != REF], rows)
            L.append("")

    Path(path).write_text("\n".join(L))
    print(f"Wrote {path}")


def write_csv(path, S, obs):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["flavour", "sample", "observable", "kind", "category", "species",
                     "pdg", "mean", "err", "n_events"])
        for fl in FLAV_ORDER:
            for s in SAMPLES:
                for o in obs:
                    i = o["i"]
                    w.writerow([FLAV_NAME[fl], s, o["key"], o["kind"], o["cat"], o["species"],
                                o["pdg"], f"{S[s]['mean'][fl, i]:.10g}",
                                f"{S[s]['err'][fl, i]:.10g}", int(S[s]["N"][fl])])
    print(f"Wrote {path}")


# ----------------------------------------------------------------------------
#  Plots
# ----------------------------------------------------------------------------
def plot_rows(rows, S, fl, title, scale, ylabel, path):
    """Publication-quality comparison plot: multiplicities + ratio to Paper.

    Plotting only: no changes to the event selection, statistics, or observables.
    Matplotlib uses the CERN/LCG LaTeX installation supplied in PATH.
    """
    import matplotlib
    matplotlib.use("Agg")

    # -------------------------------------------------------------------------
    # LaTeX / typography
    # -------------------------------------------------------------------------
    matplotlib.rcParams.update({
        "text.usetex": True,
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman"],
        "text.latex.preamble": r"\usepackage{amsmath}",
        "axes.labelsize": 11,
        "axes.titlesize": 12,
        "xtick.labelsize": 9.5,
        "ytick.labelsize": 9.5,
        "legend.fontsize": 9.5,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
        "xtick.major.size": 4,
        "ytick.major.size": 4,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "savefig.dpi": 600,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

    import matplotlib.pyplot as plt
    from matplotlib.ticker import LogLocator, NullFormatter, FuncFormatter

    idx = [o["i"] for o in rows]
    m = {s: S[s]["mean"][fl, idx] * scale for s in SAMPLES}
    e = {s: S[s]["err"][fl, idx] * scale for s in SAMPLES}

    keep = np.zeros(len(idx), dtype=bool)
    for s in SAMPLES:
        keep |= np.isfinite(m[s]) & (m[s] > 0)

    if not keep.any():
        return

    rows = [o for o, k in zip(rows, keep) if k]
    labels_raw = [o["species"] or o["key"] for o in rows]

    # Compact particle labels in LaTeX.  This keeps the tick labels readable
    # without changing the observable names stored in the analysis.
    def particle_label(name):
        special = {
            "B0": r"$B^0$",
            "B+": r"$B^+$",
            "Bs0": r"$B_s^0$",
            "Bc+": r"$B_c^+$",
            "B*0": r"$B^{*0}$",
            "B*+": r"$B^{*+}$",
            "Bs*0": r"$B_s^{*0}$",
            "Bc*+": r"$B_c^{*+}$",
            "D+": r"$D^+$",
            "D0": r"$D^0$",
            "Ds+": r"$D_s^+$",
            "D*+": r"$D^{*+}$",
            "D*0": r"$D^{*0}$",
            "Ds*+": r"$D_s^{*+}$",
            "D1(2420)+": r"$D_1(2420)^+$",
            "D1(2420)0": r"$D_1(2420)^0$",
            "D2*(2460)+": r"$D_2^*(2460)^+$",
            "D2*(2460)0": r"$D_2^*(2460)^0$",
            "Lambda_b0": r"$\Lambda_b^0$",
            "Sigma_b-": r"$\Sigma_b^-$",
            "Sigma_b0": r"$\Sigma_b^0$",
            "Sigma_b+": r"$\Sigma_b^+$",
            "Sigma_b*-": r"$\Sigma_b^{*-}$",
            "Sigma_b*0": r"$\Sigma_b^{*0}$",
            "Sigma_b*+": r"$\Sigma_b^{*+}$",
            "Xi_b-": r"$\Xi_b^-$",
            "Xi_b0": r"$\Xi_b^0$",
            "Xi_b'-": r"$\Xi_b^{\prime-}$",
            "Xi_b'0": r"$\Xi_b^{\prime0}$",
            "Xi_b*-": r"$\Xi_b^{*-}$",
            "Xi_b*0": r"$\Xi_b^{*0}$",
            "Omega_b-": r"$\Omega_b^-$",
            "Omega_b*-": r"$\Omega_b^{*-}$",
            "Lambda_c+": r"$\Lambda_c^+$",
            "Sigma_c0": r"$\Sigma_c^0$",
            "Sigma_c+": r"$\Sigma_c^+$",
            "Sigma_c++": r"$\Sigma_c^{++}$",
            "Sigma_c*0": r"$\Sigma_c^{*0}$",
            "Sigma_c*+": r"$\Sigma_c^{*+}$",
            "Sigma_c*++": r"$\Sigma_c^{*++}$",
            "Xi_c0": r"$\Xi_c^0$",
            "Xi_c+": r"$\Xi_c^+$",
            "Xi_c'0": r"$\Xi_c^{\prime0}$",
            "Xi_c'+": r"$\Xi_c^{\prime+}$",
            "Xi_c*0": r"$\Xi_c^{*0}$",
            "Xi_c*+": r"$\Xi_c^{*+}$",
            "Omega_c0": r"$\Omega_c^0$",
            "Omega_c*0": r"$\Omega_c^{*0}$",
        }
        return special.get(name, name.replace("_", r"\_"))

    labels = [particle_label(name) for name in labels_raw]
    x = np.arange(len(labels), dtype=float)

    # -------------------------------------------------------------------------
    # Figure geometry
    # -------------------------------------------------------------------------
    n = len(labels)
    width = min(11.5, max(6.8, 1.05 * n + 2.2))
    fig = plt.figure(figsize=(width, 7.0))
    gs = fig.add_gridspec(
        2, 1,
        height_ratios=(3.15, 1.45),
        hspace=0.055,
        left=0.095,
        right=0.985,
        bottom=0.20 if n > 8 else 0.16,
        top=0.93,
    )
    ax = fig.add_subplot(gs[0])
    axr = fig.add_subplot(gs[1], sharex=ax)

    # A restrained, colour-blind-friendly palette with a black reference.
    styles = {
        "Paper":   dict(marker="o", color="#222222", label="Paper"),
        "HVFL":    dict(marker="s", color="#0072B2", label="HVFL"),
        "Pythia8": dict(marker="^", color="#D55E00", label="Pythia8"),
    }

    offsets = np.linspace(-0.20, 0.20, len(SAMPLES))

    # -------------------------------------------------------------------------
    # Main multiplicity panel
    # -------------------------------------------------------------------------
    for off, s in zip(offsets, SAMPLES):
        st = styles[s]
        y = m[s][keep]
        ey = e[s][keep]
        valid = np.isfinite(y) & (y > 0)

        ax.errorbar(
            x[valid] + off, y[valid], yerr=ey[valid],
            fmt=st["marker"],
            linestyle="none",
            markersize=5.5,
            markeredgewidth=0.75,
            markeredgecolor=st["color"],
            markerfacecolor="white" if s == "Paper" else st["color"],
            ecolor=st["color"],
            elinewidth=0.9,
            capsize=2.2,
            capthick=0.8,
            label=st["label"],
            zorder=3,
        )

    ax.set_yscale("log")
    ax.set_ylabel(ylabel, labelpad=7)
    ax.set_title(title, loc="left", pad=8, fontweight="normal")

    # Log-scale grid: only major grid lines to keep the plot uncluttered.
    ax.grid(which="major", axis="y", linewidth=0.55, alpha=0.30)
    ax.grid(which="minor", axis="y", linewidth=0.35, alpha=0.12)
    ax.grid(False, axis="x")
    ax.yaxis.set_major_locator(LogLocator(base=10.0, numticks=8))
    ax.yaxis.set_minor_formatter(NullFormatter())

    ax.tick_params(axis="x", which="both", bottom=False, labelbottom=False)

    ax.legend(
        loc="upper right",
        ncol=3,
        frameon=False,
        handletextpad=0.45,
        columnspacing=1.1,
        borderaxespad=0.1,
    )

    # -------------------------------------------------------------------------
    # Ratio panel
    # -------------------------------------------------------------------------
    for off, s in zip(offsets, SAMPLES):
        if s == REF:
            continue

        st = styles[s]
        a, ea = m[s][keep], e[s][keep]
        b, eb = m[REF][keep], e[REF][keep]

        with np.errstate(divide="ignore", invalid="ignore"):
            r = np.where(b > 0, a / b, np.nan)
            er = np.where(
                (a > 0) & (b > 0),
                np.abs(r) * np.sqrt((ea / a) ** 2 + (eb / b) ** 2),
                np.nan,
            )

        valid = np.isfinite(r) & np.isfinite(er)

        axr.errorbar(
            x[valid] + off, r[valid], yerr=er[valid],
            fmt=st["marker"],
            linestyle="none",
            markersize=4.6,
            markeredgewidth=0.7,
            markeredgecolor=st["color"],
            markerfacecolor=st["color"],
            ecolor=st["color"],
            elinewidth=0.85,
            capsize=2.0,
            capthick=0.75,
            zorder=3,
        )

    axr.axhline(1.0, color="#555555", linewidth=0.8, linestyle=(0, (4, 3)), zorder=1)
    axr.set_ylabel(r"$\mathrm{ratio\ to\ Paper}$", labelpad=7)
    axr.set_xticks(x)
    axr.set_xticklabels(labels, rotation=45, ha="right", rotation_mode="anchor")

    # Keep the reference line and the useful ratio range visible while allowing
    # unusually large deviations to expand the panel rather than being clipped.
    ratio_values = []
    for s in SAMPLES:
        if s == REF:
            continue
        a, b = m[s][keep], m[REF][keep]
        with np.errstate(divide="ignore", invalid="ignore"):
            rr = np.where(b > 0, a / b, np.nan)
        ratio_values.extend(rr[np.isfinite(rr)].tolist())

    if ratio_values:
        rmax = max(1.55, float(np.nanmax(ratio_values)))
        rmin = min(0.45, float(np.nanmin(ratio_values)))
        axr.set_ylim(max(0.0, 0.85 * rmin), 1.10 * rmax)
    else:
        axr.set_ylim(0.5, 1.5)

    axr.grid(which="major", axis="y", linewidth=0.55, alpha=0.30)
    axr.grid(False, axis="x")
    axr.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))

    # Clean frame: no top/right spines, no visual box around the data.
    for a in (ax, axr):
        a.spines["top"].set_visible(False)
        a.spines["right"].set_visible(False)
        a.spines["left"].set_linewidth(0.8)
        a.spines["bottom"].set_linewidth(0.8)

    # -------------------------------------------------------------------------
    # Export: vector PDF for thesis/papers + high-resolution PNG for slides/docs
    # -------------------------------------------------------------------------
    path = Path(path)
    pdf_path = path.with_suffix(".pdf")

    fig.savefig(
        pdf_path,
        bbox_inches="tight",
        pad_inches=0.03,
        facecolor="white",
    )
    fig.savefig(
        path,
        bbox_inches="tight",
        pad_inches=0.03,
        dpi=600,
        facecolor="white",
    )
    plt.close(fig)


def make_plots(outdir, S, obs):
    try:
        import matplotlib  # noqa: F401
    except ImportError:
        print("matplotlib not available, skipping plots", file=sys.stderr)
        return

    outdir.mkdir(parents=True, exist_ok=True)

    cat_rows = [o for o in obs if o["kind"] in ("category", "weak")]

    for fl in FLAV_ORDER:
        f = FLAV_NAME[fl]

        plot_rows(
            cat_rows, S, fl,
            rf"${f}$-jets: category totals",
            1.0,
            r"$\langle n\rangle$ per event",
            outdir / f"{f}_categories.png",
        )

        for cat in CATEGORY_ORDER:
            rows = [
                o for o in obs
                if o["kind"] == "cc" and o["cat"] == cat
            ]

            plot_rows(
                rows, S, fl,
                rf"${f}$-jets: {CAT_TITLE[cat]} (+ c.c.)",
                1e3,
                r"$\langle n\rangle$ per event $\times 10^3$",
                outdir / f"{f}_{cat}.png",
            )

    print(f"Wrote plots to {outdir} (PNG + vector PDF)")

# ----------------------------------------------------------------------------
#  Main
# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", default="tables")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--max-files", type=int, default=0, help="per directory, 0 = all")
    ap.add_argument("--reprocess", action="store_true", help="ignore the cached event loop")
    ap.add_argument("--no-plots", action="store_true")
    args = ap.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    obs = build_observables()
    cache = out / "raw_counts.npz"

    loaded = None
    if cache.exists() and not args.reprocess:
        loaded = load_cache(cache, obs, args.max_files)
        if loaded:
            print(f"Using cached event-loop results from {cache} (--reprocess to redo)")
    if loaded is None:
        R, runs = run_event_loops(args, obs)
        save_cache(cache, R, runs, obs, args.max_files)
    else:
        R, runs = loaded

    S = compute(R)
    write_markdown(out / "comparison.md", S, obs, runs)
    write_csv(out / "multiplicities.csv", S, obs)
    if not args.no_plots:
        make_plots(out / "plots", S, obs)

    key = {o["key"]: o["i"] for o in obs}
    print("\nSanity check (expect about 2 per event in the matching flavour):")
    for fl, k in ((5, "weakly-decaying b hadrons"), (4, "weakly-decaying c hadrons")):
        vals = ", ".join(f"{s}: {fmt(S[s]['mean'][fl, key[k]], S[s]['err'][fl, key[k]])}"
                         for s in SAMPLES)
        print(f"  {k} in {FLAV_NAME[fl]}-events -> {vals}")


if __name__ == "__main__":
    main()
