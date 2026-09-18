"""Build the per-run list for one data-taking year from the ALEPH run database.

Inputs: the ALEPH SQLite run database (tables runs, run_detector_quality_flags)
and the per-run event-class counts of the converted data files (count_event_classes.py).
Output: run_list_<year>.csv, one row per run of the database with the detector flags,
the bookkeeping counters, the stored SICAL luminosity and the quantities needed by the
luminosity policies of src/run_list.py. The selection and the luminosity policy are
applied there, not here; this script only records the facts and the run verdict.

Verdict (column `reason`, `selected` = 1 only for `ok`):
  flags               VDET, ITC or TPC not P, or any subsystem D
  vetoed_<group>      our own veto (VETO_RUNS below), all checked on the converted files
  truncated_counters  hadronic class-16 events in the files exceed the database Z count
                      (n_z0) by more than TRUNC: the SICAL Bhabha counter, hence the stored
                      luminosity, covers only part of the run (column `truncated`, set for
                      every run whatever its flags)
  ok

Luminosity monitor: the large-angle lepton-pair rate (event class 15) per stored SICAL
nb^-1, calibrated on the PERF runs of the reference period. Fill factors are the
lepton-pair / SICAL ratio of the clean runs of each fill; `fill_rescaled` = 1 where the
factor deviates from one by more than 5 % and NSIG standard deviations.

    python build_run_list.py --db aleph.db --classes class_counts_raw_1994.json --year 1994
"""
import argparse
import collections
import csv
import json
import math
import os
import sqlite3
import sys
from datetime import date

# Our own run veto, applied on top of the database flags (1994).
# vertexing: VDET off; 26097 also a +-9.5 cm z offset per TPC drift half; 27776 sigma(d0)
#   doubled with ITC hit loss.
# tracking_quality: mean barrel TPC hit count more than 4 sigma below the neighbouring runs of
#   the same file, or the fraction of tracks with sigma(d0) < 100 um more than 4 sigma below it
#   together with an ITC/VDET hit loss; 27518, 28941, 28990 re-checked per fill: TPC hits and
#   track-class rates low, 28941 chi2/ndf ~21. 26378, 26390, 29227, 29312: raw-file comparison
#   with the runs of the same fill, ITC hits per track 6-13 % low; 26390 also sigma(d0) +30 %,
#   29312 also 49 % more events without a hadronic/lepton-pair class.
# field: solenoid off-plateau, 1.5176 T (+1.2 %) from the magnet-current readback and the track
#   curvature; compensation-coil currents at 0.4 % of their usual values.
# dedx_calibration: no dE/dx calibration banks for the run (dQdx.type 4 on every track, pads and
#   wires), so no dE/dx measurement at all; tracking normal.
VETO_RUNS = {
    1994: {
        "vertexing": (
            25521, 25522, 25527, 25528, 25530, 25531, 25907, 26097,
            26725, 26726, 26745, 26747, 27776, 29321, 29322, 29323,
        ),
        "tracking_quality": (
            25520, 25692, 25712, 25781, 25908, 25934, 25953, 26084, 26088,
            26378, 26390, 26746, 26854, 27518, 28067, 28206, 28530, 28941,
            28990, 29116, 29227, 29312, 29539, 29545, 29601, 29977, 30079,
            30080, 30086, 30095, 30159, 30160, 30177, 30188, 30322, 30340,
            30347,
        ),
        "field": (28125,),
        "dedx_calibration": (26420, 27750, 29667),
    },
}
REFERENCE_PERIOD = {1994: ("1994-07-01", "1994-09-30")}
SUBSYSTEMS = ("vdet", "itc", "tpc", "ecal", "hcal", "sical", "lcal", "trigger", "daq", "run_coord")
TRACKING = ("vdet", "itc", "tpc")
LEPTON_PAIR_CLASS = 15
HADRONIC_CLASS = 16

COLUMNS = ["run", "fill", "date", "ecm_gev", "run_quality", *SUBSYSTEMS,
           "n_events", "n_z0", "n_bhabha", "lumi_sical_nb", "lumi_status",
           "class15", "class16", "truncated", "reason", "selected",
           "fill_factor", "fill_factor_err", "fill_rescaled", "lumi_leptonpair_nb"]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", required=True, help="ALEPH run database (SQLite)")
    ap.add_argument("--classes", required=True, help="per-run event-class counts JSON from count_event_classes.py")
    ap.add_argument("--year", type=int, default=1994)
    ap.add_argument("--trunc", type=float, default=1.1, help="class-16 / n_z0 above which the counters are truncated")
    ap.add_argument("--nsig", type=float, default=4.0, help="significance for a fill rescaling")
    ap.add_argument("--out", default=None, help="output CSV (default: run_list_<year>.csv next to this script)")
    a = ap.parse_args()
    out = a.out or os.path.join(os.path.dirname(os.path.abspath(__file__)), f"run_list_{a.year}.csv")

    veto = {r: g for g, rs in VETO_RUNS.get(a.year, {}).items() for r in rs}
    classes = json.load(open(a.classes))["per_run"]

    def cls(run, k):
        return classes.get(str(run), {}).get("cls", {}).get(str(k), 0)

    con = sqlite3.connect(f"file:{a.db}?mode=ro", uri=True)
    con.row_factory = sqlite3.Row
    rows = {r["run_number"]: r for r in con.execute(
        "select r.*, f.* from runs r join run_detector_quality_flags f using(run_number) where r.year = ?", (a.year,))}
    if not rows:
        sys.exit(f"no runs of {a.year} in {a.db}")

    def flags_ok(x):
        return all(x[f"{s}_qual"] == "P" for s in TRACKING) and all(x[f"{s}_qual"] != "D" for s in SUBSYSTEMS)

    def lumi_known(x):
        return x["lumi_status"] == "OK" and (x["lumi_nb"] or 0) > 0

    # lepton-pair rate per nb^-1 on the PERF runs of the reference period
    if a.year not in REFERENCE_PERIOD:
        sys.exit(f"no lepton-pair reference period defined for {a.year} (REFERENCE_PERIOD)")
    d0, d1 = REFERENCE_PERIOD[a.year]
    # truncated bookkeeping counters (checked for every run: it matters for the reference rate and the fill factors too)
    truncated = {r for r, x in rows.items()
                 if lumi_known(x) and cls(r, HADRONIC_CLASS) >= 50 and cls(r, HADRONIC_CLASS) > a.trunc * max(x["n_z0"] or 0, 1)}
    ref = [r for r, x in rows.items() if r not in veto and r not in truncated and x["run_quality"] == "PERF" and lumi_known(x) and x["run_date"] and d0 <= x["run_date"] <= d1]
    L_ref = sum(rows[r]["lumi_nb"] for r in ref)
    rate15 = sum(cls(r, LEPTON_PAIR_CLASS) for r in ref) / L_ref

    # verdicts
    reason = {}
    for r, x in rows.items():
        if not flags_ok(x):
            reason[r] = "flags"
        elif r in veto:
            reason[r] = "vetoed_" + veto[r]
        elif r in truncated:
            reason[r] = "truncated_counters"
        else:
            reason[r] = "ok"

    # fill factors from the lumi-known runs of each fill, whatever their detector flags (the
    # monitor tests the SICAL bookkeeping), without the vetoed and truncated runs
    fills = collections.defaultdict(lambda: {"L": 0.0, "c15": 0, "nb": 0})
    for r, x in rows.items():
        if r in veto or r in truncated or not lumi_known(x):
            continue
        f = fills[x["fill_number"]]
        f["L"] += x["lumi_nb"]
        f["c15"] += cls(r, LEPTON_PAIR_CLASS)
        f["nb"] += x["n_bhabha"] or 0
    factor = {}
    for fn, f in fills.items():
        if f["L"] > 0 and f["c15"] > 0:
            ratio = f["c15"] / rate15 / f["L"]
            err = ratio * math.sqrt(1 / f["c15"] + 1 / max(f["nb"], 1))
            factor[fn] = (ratio, err, abs(ratio - 1) > 0.05 and abs(ratio - 1) > a.nsig * err)

    table = []
    for r in sorted(rows):
        x = rows[r]
        fr, fe, resc = factor.get(x["fill_number"], (1.0, 0.0, False))
        table.append({
            "run": r, "fill": x["fill_number"], "date": x["run_date"], "ecm_gev": x["ecm_gev"], "run_quality": x["run_quality"],
            **{s: x[f"{s}_qual"] for s in SUBSYSTEMS},
            "n_events": x["n_events"], "n_z0": x["n_z0"], "n_bhabha": x["n_bhabha"],
            "lumi_sical_nb": x["lumi_nb"] if lumi_known(x) else "", "lumi_status": x["lumi_status"],
            "class15": cls(r, LEPTON_PAIR_CLASS), "class16": cls(r, HADRONIC_CLASS),
            "truncated": int(r in truncated), "reason": reason[r], "selected": int(reason[r] == "ok"),
            "fill_factor": f"{fr:.4f}", "fill_factor_err": f"{fe:.4f}", "fill_rescaled": int(resc),
            "lumi_leptonpair_nb": f"{cls(r, LEPTON_PAIR_CLASS) / rate15:.3f}",
        })

    header = [
        f"# ALEPH {a.year} run list, built {date.today().isoformat()} by build_run_list.py",
        f"# database: {os.path.basename(a.db)} (ALEPH run database, SQLite); event-class counts: {os.path.basename(a.classes)} (count_event_classes.py)",
        f"# selection: VDET, ITC, TPC = P and no subsystem D; our veto {len(veto)} runs; truncated counters: class16 / n_z0 > {a.trunc}",
        f"# lepton-pair monitor (class {LEPTON_PAIR_CLASS}): {rate15:.4f} per nb^-1 on {len(ref)} PERF runs {d0}..{d1}, {L_ref / 1e3:.2f} pb^-1 SICAL",
        f"# fill factor = lepton-pair / SICAL luminosity of the fill's clean runs; fill_rescaled = 1 when |factor - 1| > 0.05 and > {a.nsig} sigma",
        "# luminosities in nb^-1; the luminosity policy is applied by src/run_list.py",
    ]
    with open(out, "w", newline="") as fh:
        fh.write("\n".join(header) + "\n")
        w = csv.DictWriter(fh, fieldnames=COLUMNS, lineterminator="\n")
        w.writeheader()
        w.writerows(table)

    n = collections.Counter(t["reason"] for t in table)
    print(f"{len(table)} runs of {a.year} -> {out}")
    for k in sorted(n):
        print(f"  {k:<26} {n[k]:>5}")
    print("rescaled fills:", {fn: round(v[0], 3) for fn, v in sorted(factor.items()) if v[2]})


if __name__ == "__main__":
    main()
