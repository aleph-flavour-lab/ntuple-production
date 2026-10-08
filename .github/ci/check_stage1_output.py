#!/usr/bin/env python3
"""Checks one stage1 output of the CI and prints what the synthetic input exercised.

Fails (exit status 1) when the file has no `events` tree, no entries or no
`eventsProcessed`. The other numbers are printed for information only: they show
which parts of stage1 found something to do (vertex fits, finders, dE/dx) and are
not compared with any reference.
Usage: check_stage1_output.py OUTPUT.root [--branches-out FILE]
"""
import argparse
import sys

import awkward as ak
import numpy as np
import uproot


def _count(t, name):
    return int(ak.sum(t[name].array())) if name in t else "n/a"


def _coverage(t):
    """(label, value) rows; a branch missing from the output gives 'n/a'."""
    rows = []

    def row(label, fn):
        try:
            rows.append((label, fn()))
        except KeyError:
            rows.append((label, "n/a"))

    n = t.num_entries

    def distinct(name):
        return len(np.unique(ak.to_numpy(ak.flatten(t[name].array(), axis=None))))
    row("runs / beam-spot centres", lambda: f"{distinct('run_number')} / {distinct('Beamspot_x')}")
    row("events with a usable PV (pv_good)", lambda: f"{int(ak.sum(t['pv_good'].array()))} / {n}")
    row("PV fit failed / pruning failed / trivial (events)",
        lambda: " / ".join(str(int(ak.sum(t[b].array() == v)))
                           for b, v in (("pv_converged", 0), ("pv_split_converged", 0), ("pv_trivial", 1))))
    row("tracks per event: all / baseline / PV fit / primary",
        lambda: " / ".join(f"{ak.mean(t[b].array()):.1f}" for b in
                           ("n_tracks_all", "n_tracks_sel", "n_tracks_sel_vertexfit", "n_primary_tracks")))
    row("legacy SV / V0 candidates", lambda: f"{_count(t, 'n_sv_event')} / {_count(t, 'n_v0_event')}")

    def v0n():
        tight = t["v0n_tight"].array()
        pdg = t["v0n_pdg"].array()
        mass = t["v0n_invM"].array()
        ks = mass[(tight == 1) & (pdg == 310)]
        lam = mass[(tight == 1) & (pdg == 3122)]
        out = f"{int(ak.sum(ak.num(pdg)))} ({int(ak.sum(tight))} tight)"
        if ak.sum(ak.num(ks)):
            out += f", tight K0s mean mass {ak.mean(ks):.4f} GeV"
        if ak.sum(ak.num(lam)):
            out += f", tight Lambda {int(ak.sum(ak.num(lam)))}"
        return out
    row("V0 module candidates (v0n)", v0n)
    row("SV module candidates (svn)", lambda: str(_count(t, "n_svn_event")))
    row("phi->KK candidates", lambda: str(_count(t, "n_phikk_event")))
    row("D* candidates / D0 fits", lambda: f"{_count(t, 'n_dstar_event')} / {_count(t, 'n_d0fits_event')}")

    def dedx():
        q = ak.flatten(t["pfcand_charge"].array(), axis=None)
        out = []
        for leg in ("pads", "wires"):
            v = ak.flatten(t[f"pfcand_dEdx_{leg}_value"].array(), axis=None)
            ch = q != 0
            out.append(f"{leg} {float(ak.sum((v > 0) & ch)) / max(float(ak.sum(ch)), 1.0):.2f}")
        return ", ".join(out)
    row("charged constituents with a valid dE/dx", dedx)
    row("jet constituents per event", lambda: f"{ak.mean(ak.sum(t['jet_nconst'].array(), axis=1)):.1f}")
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("output")
    ap.add_argument("--branches-out", help="write the sorted branch list here")
    args = ap.parse_args()

    problems = []
    try:
        f = uproot.open(args.output)
    except Exception as exc:  # noqa: BLE001 - any failure to open is a failed check
        print(f"FAIL cannot open {args.output}: {exc}")
        return 1
    if "events" not in f:
        print(f"FAIL {args.output} has no events tree")
        return 1
    t = f["events"]
    n = t.num_entries
    branches = sorted(t.keys())
    if n == 0:
        problems.append("the events tree has no entries")
    if "eventsProcessed" not in f:
        problems.append("no eventsProcessed")
    processed = f["eventsProcessed"].member("fVal") if "eventsProcessed" in f else "n/a"
    if args.branches_out:
        with open(args.branches_out, "w") as fh:
            fh.write("\n".join(branches) + "\n")

    print(f"{args.output}: {n} entries ({processed} events processed), {len(branches)} branches")
    if n:
        for label, value in _coverage(t):
            print(f"  {label:52s} {value}")
    for p in problems:
        print(f"FAIL {p}")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
