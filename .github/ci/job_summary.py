#!/usr/bin/env python3
"""Markdown summary of a stage1 CI run, for the GitHub job summary page.

Reads what .github/ci/run_stage1_ci.sh leaves in its work directories: the step status
lines (WORK/logs/<step>.status), the output checks (WORK/output/check_<mode>.txt) and,
for a comparison with a base (main for a pull request), the comparison results
(WORK/output/compare_<mode>.json, written by compare_stage1_outputs.py).
The checkout under test is called "this PR" in a pull-request run (GITHUB_EVENT_NAME, set
by GitHub Actions), "this run" otherwise.
Python standard library only, so that it runs outside the Key4hep environment.
Usage: job_summary.py WORK [--base-work BASE_WORK] [--base-label LABEL]
"""
import argparse
import json
import os

STEPS = ("build", "inputs", "data", "mc", "compare")
MODES = (("data", "data"), ("mc", "MC"))
OPEN_TABLE_ROWS = 25   # longer tables of changed branches start collapsed
THIS = "this PR" if os.environ.get("GITHUB_EVENT_NAME") == "pull_request" else "this run"


def read(path):
    try:
        with open(path) as fh:
            return fh.read()
    except OSError:
        return None


def step_status(work, step):
    """'PASS (137 s)' from a status line 'data    PASS   137s', or None if the step did not run."""
    line = read(os.path.join(work, "logs", f"{step}.status"))
    if not line:
        return None
    parts = line.split()
    return f"{parts[1]} ({parts[2].rstrip('s')} s)" if len(parts) >= 3 else line.strip()


def fmt(x):
    return "-" if x is None else f"{x:.3g}"


def status_table(work, base_work, base_label):
    out = []
    if base_work:
        out += [f"| step | {THIS} | {base_label} |", "|---|---|---|"]
    else:
        out += ["| step | status |", "|---|---|"]
    for step in STEPS:
        head = step_status(work, step)
        if base_work:
            if step in ("build", "inputs", "compare"):
                base = step_status(base_work, step) or ("" if step == "compare" else "shared")
            else:
                base = step_status(base_work, step) or "not run"
            if head is None and not base:
                continue
            out.append(f"| {step} | {head or 'not run'} | {base} |")
        elif head is not None:
            out.append(f"| {step} | {head} |")
    return out


def names(lst):
    return ", ".join(f"`{n}`" for n in lst) if lst else "none"


def comparison(work, base_label):
    results, reasons = {}, {}
    for mode, _ in MODES:
        js = read(os.path.join(work, "output", f"compare_{mode}.json"))
        if js:
            results[mode] = json.loads(js)
        else:
            txt = read(os.path.join(work, "output", f"compare_{mode}.txt")) or "not compared"
            reasons[mode] = txt.strip().splitlines()[0] if txt.strip() else "not compared"
    out = [f"### Output changes with respect to {base_label}, on the synthetic input", ""]
    out.append("Informational: differences are listed here and do not fail the job.")
    tol = next(iter(results.values()))["tolerance"] if results else None
    if tol and (tol["rtol"] or tol["atol"]):
        out.append(f"Values compared with rtol = {tol['rtol']}, atol = {tol['atol']}.")
    elif tol:
        out.append("Values compared exactly, event by event (events matched by run and event number).")
    out += ["", "| | " + " | ".join(label for _, label in MODES) + " |", "|---|" + "---|" * len(MODES)]

    def row(label, fn, first=False):   # a mode without comparison: its reason in the first row
        cells = []
        for mode, _ in MODES:
            cells.append(fn(results[mode]) if mode in results else reasons[mode] if first else "-")
        out.append(f"| {label} | " + " | ".join(cells) + " |")

    row("result", lambda r: "**identical**" if r["identical"] else "**different**", first=True)
    row("events compared", lambda r: str(r["events"]["matched"]))
    row(f"events only in {base_label} / only in {THIS}",
        lambda r: f"{r['events']['only_ref']} / {r['events']['only_new']}")
    row("branches added", lambda r: str(len(r["added"])))
    row("branches removed", lambda r: str(len(r["removed"])))
    row("branches changed", lambda r: f"{len(r['changed'])} of {r['branches']['common']}")
    out.append("")

    for key, title in (("added", "Branches added"), ("removed", "Branches removed")):
        lists = {mode: res[key] for mode, res in results.items() if res[key]}
        if not lists:
            continue
        if len(lists) == len(results) and len({tuple(v) for v in lists.values()}) == 1:
            out.append(f"**{title}** ({' and '.join(l for m, l in MODES if m in results)}): "
                       f"{names(next(iter(lists.values())))}")
        else:
            for mode, label in MODES:
                if mode in lists:
                    out.append(f"**{title}** ({label}): {names(lists[mode])}")
        out.append("")

    for mode, label in MODES:
        changed = results.get(mode, {}).get("changed")
        if not changed:
            continue
        groups = {}
        for c in changed:
            prefix = c["branch"].split("_")[0] + "_*" if "_" in c["branch"] else c["branch"]
            groups[prefix] = groups.get(prefix, 0) + 1
        out.append(f"**Branches changed ({label})**, by name: " + ", ".join(
            f"`{g}` {n}" for g, n in sorted(groups.items(), key=lambda x: (-x[1], x[0]))))
        out.append("")
        opened = " open" if len(changed) <= OPEN_TABLE_ROWS else ""
        out.append(f"<details{opened}><summary>All {len(changed)} changed branches ({label})</summary>")
        out += ["", "| branch | events changed | fraction | largest abs. difference | largest rel. difference | note |",
                "|---|---|---|---|---|---|"]
        for c in changed:
            out.append(f"| `{c['branch']}` | {c['events']} | {c['fraction']:.3f} | {fmt(c['max_abs'])} | "
                       f"{fmt(c['max_rel'])} | {c['note']} |")
        out += ["", "</details>", ""]
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("work", help="work directory of the checkout under test")
    ap.add_argument("--base-work", help="work directory of the base, when there was a comparison")
    ap.add_argument("--base-label", default="main", help="name of the base in the tables")
    args = ap.parse_args()
    base_work = args.base_work if args.base_work and os.path.isdir(args.base_work) else None

    out = ["### stage1 on synthetic input", ""]
    first = (read(os.path.join(args.work, "summary.txt")) or "").splitlines()
    if first:
        out += [first[0], ""]
    out += status_table(args.work, base_work, args.base_label) + [""]
    if base_work:
        out += comparison(args.work, args.base_label)
    checks = [read(os.path.join(args.work, "output", f"check_{m}.txt")) for m, _ in MODES]
    if any(checks):
        out += [f"<details><summary>What the synthetic input exercised ({THIS})</summary>", "", "```"]
        out += [c.rstrip() for c in checks if c]
        out += ["```", "", "</details>"]
    print("\n".join(out))


if __name__ == "__main__":
    main()
