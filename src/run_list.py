"""Run selection and luminosity of the ALEPH data from data/lumi/run_list_<year>.csv
(built by data/lumi/build_run_list.py from the ALEPH run database).

stage1.py keeps only the selected runs; the plotting configuration takes the luminosity
from luminosity_pb(). Override the file with $ALEPH_RUN_LIST_<year>.

Luminosity policies (per selected run, nb^-1 in the file, pb^-1 here):
  nominal   stored SICAL luminosity; in fills whose SICAL bookkeeping is off (fill_rescaled = 1)
            multiplied by the fill's lepton-pair factor; runs without a stored luminosity
            valued by their large-angle lepton-pair rate
  sical     stored SICAL luminosity only, unscaled; runs without one are dropped
  sical_clean_fills   as sical, and the runs of the rescaled fills are dropped too
"""
import csv
import os

DEFAULT_POLICY = "nominal"
POLICIES = ("nominal", "sical", "sical_clean_fills")
_cache = {}


def run_list_file(year):
    return os.environ.get(f"ALEPH_RUN_LIST_{year}",
                          os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "lumi", f"run_list_{year}.csv"))


def load(year):
    """All rows of the year's list, keyed by run number (empty dict if there is no list)."""
    year = int(year)
    if year not in _cache:
        path = run_list_file(year)
        rows = {}
        if os.path.exists(path):
            with open(path) as fh:
                rows = {int(r["run"]): r for r in csv.DictReader(line for line in fh if not line.startswith("#"))}
        _cache[year] = rows
    return _cache[year]


def run_luminosities(year, policy=DEFAULT_POLICY):
    """{run: luminosity in pb^-1} of the runs kept under the policy."""
    if policy not in POLICIES:
        raise ValueError(f"unknown luminosity policy {policy!r}; choose from {POLICIES}")
    rows = load(year)
    if not rows:
        raise FileNotFoundError(f"no run list for {year}: {run_list_file(year)}")
    out = {}
    for run, r in rows.items():
        if r["selected"] != "1":
            continue
        sical = float(r["lumi_sical_nb"]) if r["lumi_sical_nb"] else None
        rescaled = r["fill_rescaled"] == "1"
        if policy == "nominal":
            lumi = float(r["lumi_leptonpair_nb"]) if sical is None else sical * (float(r["fill_factor"]) if rescaled else 1.0)
        elif sical is None or (policy == "sical_clean_fills" and rescaled):
            continue
        else:
            lumi = sical
        out[run] = lumi / 1e3
    return out


def has_list(year):
    return bool(load(year))


def good_runs(year, policy=DEFAULT_POLICY):
    """Run numbers kept under the policy."""
    return frozenset(run_luminosities(year, policy))


def luminosity_pb(year, policy=DEFAULT_POLICY):
    return sum(run_luminosities(year, policy).values())


def summary(year, policy=DEFAULT_POLICY, exclude=()):
    """One line: kept runs (minus `exclude`), luminosity, and its composition."""
    rows = load(year)
    lumis = {run: lumi for run, lumi in run_luminosities(year, policy).items() if run not in set(exclude)}
    parts = {"SICAL": [0, 0.0], "SICAL x fill factor": [0, 0.0], "lepton-pair rate": [0, 0.0]}
    for run, lumi in lumis.items():
        r = rows[run]
        key = "lepton-pair rate" if not r["lumi_sical_nb"] else ("SICAL x fill factor" if r["fill_rescaled"] == "1" and policy == "nominal" else "SICAL")
        parts[key][0] += 1
        parts[key][1] += lumi
    detail = ", ".join(f"{k} {n} runs {L:.2f}" for k, (n, L) in parts.items() if n)
    return f"run list {year} policy {policy}: {len(lumis)} of {len(rows)} runs, {sum(lumis.values()):.2f} pb^-1 ({detail})"


if __name__ == "__main__":
    import sys
    year = sys.argv[1] if len(sys.argv) > 1 else 1994
    for policy in POLICIES:
        print(summary(year, policy))
