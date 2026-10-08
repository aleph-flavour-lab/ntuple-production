#!/usr/bin/env python3
"""Compares two stage1 outputs event by event, matched by (run_number, event_number).

REF is the reference output (in the CI: main), NEW the output under test (the pull
request). The multithreaded event loop writes the events in a different order from run
to run, so the comparison never relies on the entry order. Reported:
  - branches only in NEW (added) and only in REF (removed);
  - events only in one of the two outputs (an event selection that changed);
  - for every common branch, the matched events whose values differ, with the fraction
    of matched events and the largest absolute and relative difference |new - ref| and
    |new - ref| / max(|new|, |ref|) over the values that can be paired. Events whose
    list branches have a different number of entries count as differing.
Values are compared exactly by default (NaN equals NaN); with --rtol/--atol a value
differs when |new - ref| > atol + rtol * |ref|.
Exit status: 0 identical, 1 differences found, 2 the comparison could not be made.
Usage: compare_stage1_outputs.py REF.root NEW.root [--rtol R] [--atol A] [--json FILE] [--jobs N]
"""
import argparse
import json
import multiprocessing
import os
import sys

import awkward as ak
import numpy as np
import uproot


def _leaf_equal(new, ref, rtol, atol):
    new = ak.to_numpy(new)
    ref = ak.to_numpy(ref)
    if new.dtype.kind == "f" or ref.dtype.kind == "f":
        new = new.astype(np.float64)
        ref = ref.astype(np.float64)
        if rtol == 0 and atol == 0:
            return (new == ref) | (np.isnan(new) & np.isnan(ref))
        return np.isclose(new, ref, rtol=rtol, atol=atol, equal_nan=True)
    return new == ref


def _per_event(fn, new, ref):
    """Applies a per-value test fn at any list depth; an event passes when its lists have
    the same lengths and every value passes (numpy bool array, one entry per event)."""
    if new.ndim == 1:
        return fn(new, ref)
    ok = ak.to_numpy(ak.num(new, axis=1) == ak.num(ref, axis=1))
    idx = np.nonzero(ok)[0]
    if len(idx):
        sn, sr = new[idx], ref[idx]
        inner = _per_event(fn, ak.flatten(sn, axis=1), ak.flatten(sr, axis=1))
        ok[idx] = ak.to_numpy(ak.all(ak.unflatten(inner, ak.num(sn, axis=1)), axis=1))
    return ok


def event_equal(new, ref, rtol=0.0, atol=0.0):
    """Per-event equality of two equally long arrays of any list depth."""
    return _per_event(lambda n, r: _leaf_equal(n, r, rtol, atol), new, ref)


def largest_differences(new, ref):
    """(max |new-ref|, max |new-ref|/max(|new|,|ref|), number of NaN/inf mismatches) over
    the events whose lists have the same lengths in both outputs."""
    idx = np.nonzero(_same_shape(new, ref))[0]
    if not len(idx):
        return None, None, 0
    fn = ak.to_numpy(ak.flatten(new[idx], axis=None)).astype(np.float64)
    fr = ak.to_numpy(ak.flatten(ref[idx], axis=None)).astype(np.float64)
    finite = np.isfinite(fn) & np.isfinite(fr)
    both_nan = np.isnan(fn) & np.isnan(fr)
    nonfinite_mismatch = int(np.sum(~finite & ~both_nan & (fn != fr)))
    d = np.abs(fn[finite] - fr[finite])
    if not len(d):
        return None, None, nonfinite_mismatch
    scale = np.maximum(np.abs(fn[finite]), np.abs(fr[finite]))
    rel = np.divide(d, scale, out=np.zeros_like(d), where=scale > 0)
    return float(d.max()), float(rel.max()), nonfinite_mismatch


def _scalar_column(t, name):
    """A per-event number; EventHeader-derived branches hold it as a one-element vector."""
    arr = t[name].array()
    return ak.to_numpy(ak.firsts(arr) if arr.ndim > 1 else arr)


def _keys(t):
    run = _scalar_column(t, "run_number")
    evt = _scalar_column(t, "event_number")
    return [(int(r), int(e)) for r, e in zip(run, evt)]


_W = {}   # per worker process: the two trees and the matched entries


def _init_worker(ref_path, new_path, i_ref, i_new, rtol, atol):
    _W.update(tr=uproot.open(ref_path)["events"], tn=uproot.open(new_path)["events"],
              i_ref=i_ref, i_new=i_new, rtol=rtol, atol=atol)


def _same_shape(new, ref):
    return _per_event(lambda n, r: np.ones(len(n), dtype=bool), new, ref)


def _compare_branch(name):
    """One common branch over the matched events: a row of the result (events = 0: unchanged)."""
    n_match = len(_W["i_new"])
    row = {"branch": name, "events": 0, "fraction": 0.0, "max_abs": None, "max_rel": None, "note": ""}
    try:
        new = _W["tn"][name].array()[_W["i_new"]]
        ref = _W["tr"][name].array()[_W["i_ref"]]
        if new.ndim != ref.ndim:
            row.update(events=n_match, note="list depth differs")
        else:
            nbad = int(np.sum(~event_equal(new, ref, _W["rtol"], _W["atol"])))
            if nbad:
                row["events"] = nbad
                nshape = int(np.sum(~_same_shape(new, ref)))
                notes = []
                if nshape:
                    notes.append(f"{nshape} events with a different number of entries")
                if nbad > nshape:   # some events differ in values, not only in length
                    try:
                        row["max_abs"], row["max_rel"], nonfinite = largest_differences(new, ref)
                    except Exception as exc:  # noqa: BLE001 - the counts above stay valid
                        nonfinite = 0
                        notes.append(f"largest difference not computed: {type(exc).__name__}")
                    if nonfinite:
                        notes.append(f"{nonfinite} NaN/inf values differ")
                row["note"] = "; ".join(notes)
    except Exception as exc:  # noqa: BLE001 - reported per branch, the rest is still compared
        row.update(events=n_match, note=f"not comparable: {type(exc).__name__}: {' '.join(str(exc).split())[:200]}")
    if row["events"]:
        row["fraction"] = row["events"] / n_match if n_match else 0.0
    return row


def compare(ref_path, new_path, rtol=0.0, atol=0.0, jobs=1):
    """The comparison as a dict (written by --json; the CI job summary is made from it).
    The branches are read and compared in `jobs` processes."""
    with uproot.open(ref_path) as fr, uproot.open(new_path) as fn:
        tr, tn = fr["events"], fn["events"]
        kr, kn = _keys(tr), _keys(tn)
        br, bn = set(tr.keys()), set(tn.keys())
    for label, keys in (("REF", kr), ("NEW", kn)):
        if len(set(keys)) != len(keys):
            raise ValueError(f"(run_number, event_number) is not unique in {label}")
    pos_r = {k: i for i, k in enumerate(kr)}
    matched = [(i, pos_r[k]) for i, k in enumerate(kn) if k in pos_r]
    i_new = np.array([m[0] for m in matched], dtype=np.int64)
    i_ref = np.array([m[1] for m in matched], dtype=np.int64)
    n_match = len(matched)

    names = sorted(br & bn)
    init = (ref_path, new_path, i_ref, i_new, rtol, atol)
    if jobs > 1:
        with multiprocessing.get_context("fork").Pool(jobs, initializer=_init_worker, initargs=init) as pool:
            rows = pool.map(_compare_branch, names, chunksize=4)
    else:
        _init_worker(*init)
        rows = [_compare_branch(n) for n in names]
    changed = [r for r in rows if r["events"]]

    result = {
        "ref": ref_path, "new": new_path,
        "tolerance": {"rtol": rtol, "atol": atol},
        "events": {"ref": len(kr), "new": len(kn), "matched": n_match,
                   "only_ref": len(kr) - n_match, "only_new": len(kn) - n_match},
        "branches": {"ref": len(br), "new": len(bn), "common": len(br & bn)},
        "added": sorted(bn - br), "removed": sorted(br - bn), "changed": changed,
    }
    e = result["events"]
    result["identical"] = not (changed or result["added"] or result["removed"] or e["only_ref"] or e["only_new"])
    return result


def _fmt(x):
    return "-" if x is None else f"{x:.3g}"


def print_text(r):
    e, b = r["events"], r["branches"]
    tol = r["tolerance"]
    print(f"REF: {r['ref']} ({e['ref']} events, {b['ref']} branches)")
    print(f"NEW: {r['new']} ({e['new']} events, {b['new']} branches)")
    print(f"matched events: {e['matched']}; only in REF: {e['only_ref']}; only in NEW: {e['only_new']}")
    print(f"branches added (only in NEW):   {r['added'] or 'none'}")
    print(f"branches removed (only in REF): {r['removed'] or 'none'}")
    how = "exact" if tol["rtol"] == 0 and tol["atol"] == 0 else f"rtol={tol['rtol']} atol={tol['atol']}"
    print(f"common branches changed ({how}): {len(r['changed'])} of {b['common']}")
    if r["changed"]:
        print(f"  {'branch':45s} {'events':>7s} {'fraction':>8s} {'max abs':>10s} {'max rel':>10s}")
    for c in r["changed"]:
        print(f"  {c['branch']:45s} {c['events']:7d} {c['fraction']:8.3f} {_fmt(c['max_abs']):>10s} "
              f"{_fmt(c['max_rel']):>10s}  {c['note']}")
    print("IDENTICAL" if r["identical"] else "DIFFERENT")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("ref", help="reference stage1 output (e.g. main)")
    ap.add_argument("new", help="stage1 output under test (e.g. the pull request)")
    ap.add_argument("--rtol", type=float, default=0.0, help="relative tolerance (default 0: exact)")
    ap.add_argument("--atol", type=float, default=0.0, help="absolute tolerance (default 0: exact)")
    ap.add_argument("--json", help="also write the result to this JSON file")
    ap.add_argument("--jobs", type=int, default=len(os.sched_getaffinity(0)),
                    help="processes reading and comparing the branches (default: the available cores)")
    args = ap.parse_args()
    try:
        result = compare(args.ref, args.new, args.rtol, args.atol, max(1, args.jobs))
    except Exception as exc:  # noqa: BLE001 - any failure here means no comparison
        print(f"comparison not possible: {type(exc).__name__}: {' '.join(str(exc).split())[:300]}")
        return 2
    print_text(result)
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(result, fh, indent=1)
    return 0 if result["identical"] else 1


if __name__ == "__main__":
    sys.exit(main())
