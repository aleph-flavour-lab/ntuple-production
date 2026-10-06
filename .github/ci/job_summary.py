#!/usr/bin/env python3
"""Markdown summary of a stage1 CI run, for the GitHub job summary page.

Reads what .github/ci/run_stage1_ci.sh leaves in its work directories: the step status
lines (WORK/logs/<step>.status), the output checks (WORK/output/check_<mode>.txt) and,
for a comparison with a base (main for a pull request), the comparison results
(WORK/output/compare_<mode>.json, written by compare_stage1_outputs.py). For a failed step,
up to five error lines (always the last one) of the file that tells why are quoted: its log
(WORK/logs/<step>.log, stage1_<mode>.log), the output check or the comparison output.
The checkout under test is called "this PR" in a pull-request run (GITHUB_EVENT_NAME, set
by GitHub Actions), "this run" otherwise.
Python standard library only, so that it runs outside the Key4hep environment.
Usage: job_summary.py WORK [--base-work BASE_WORK] [--base-label LABEL]
"""
import argparse
import json
import os
import re
import unicodedata

STEPS = ("build", "inputs", "data", "mc", "compare")
MODES = (("data", "data"), ("mc", "MC"))
OPEN_TABLE_ROWS = 25   # longer tables of changed branches start collapsed
THIS = "this PR" if os.environ.get("GITHUB_EVENT_NAME") == "pull_request" else "this run"
MAX_ERROR_LINES = 5    # error lines quoted for a failed step
MAX_LINE = 200         # characters per quoted line
# lines of a log that tell why a step failed: diagnostics of the compiler and of the
# interpreter that compiles the analysis (cling), errors of FCCAnalyses, ROOT and CMake, the
# last line of a Python traceback, crashes, the FAIL lines of check_stage1_output.py, a
# comparison that could not be made, and the ERROR lines written by run_stage1_ci.sh
ERROR_LINE = re.compile(r"(?i:\berror:)|\bERROR\b|^(Error|Fatal|SysError) in <|^CMake Error"
                        r"|^[A-Za-z_][\w.]*(Error|Exception)(: |$)|^what\(\):|^terminate called"
                        r"|\*\*\* Break \*\*\*|(?i:segmentation (fault|violation))|^FAIL |^comparison not possible")
ESCAPE = re.compile(r"\x1b(\[[0-9;?]*[ -/]*[@-~]|[@-Z\\-_])")   # terminal colours and the like
UP = re.compile(r"/(?!\.\.?/)[^/\s'\"]+/\.\.(?=/)")            # 'dir/..' in a path


def read(path):
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
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


def plain(line):
    """A log line without terminal escapes, other control or format characters and outer spaces."""
    line = ESCAPE.sub("", line).replace("\t", " ")
    return "".join(c for c in line if unicodedata.category(c) not in ("Cc", "Cf")).strip()


def shorten(line):
    """'a/b/../c' as 'a/c' (the analysis headers are included as .github/ci/../../src/...), and
    at most MAX_LINE characters."""
    prev = None
    while prev != line:
        prev, line = line, UP.sub("", line)
    return line if len(line) <= MAX_LINE else line[:MAX_LINE - 3] + "..."


def error_lines(text):
    """The distinct error lines of a log (lines differing only in numbers count once), at most
    MAX_ERROR_LINES: the first ones, and always the last one (run_stage1_ci.sh appends its own
    reason, such as a timeout, at the end of the log). A line ending with ':' gets the indented
    line right after it, which holds the message (as in FCCAnalyses' 'ERROR: During the
    execution of the analysis file exception occurred:')."""
    lines = text.splitlines()
    found, last_key, last_line = {}, None, None   # key -> first line with that key, in log order
    for i, raw in enumerate(lines):
        line = ESCAPE.sub("", raw).strip()
        if not ERROR_LINE.search(line):
            continue
        line = plain(line)
        nxt = ESCAPE.sub("", lines[i + 1]) if i + 1 < len(lines) else ""
        if line.endswith(":") and nxt[:1].isspace() and not ERROR_LINE.search(nxt.strip()):
            line += " " + plain(nxt)
        last_key, last_line = re.sub(r"[0-9]+", "0", line), line
        found.setdefault(last_key, line)
    keys = list(found)
    if len(keys) <= MAX_ERROR_LINES or last_key in keys[:MAX_ERROR_LINES - 1]:
        return [shorten(found[k]) for k in keys[:MAX_ERROR_LINES]]
    return [shorten(found[k]) for k in keys[:MAX_ERROR_LINES - 1]] + [shorten(last_line)]


def last_lines(text):
    """The last MAX_ERROR_LINES lines of a log that are not blank."""
    out = []
    for raw in reversed(text.splitlines()):
        line = plain(raw)
        if line:
            out.insert(0, shorten(line))
            if len(out) == MAX_ERROR_LINES:
                break
    return out


def code_block(lines):
    """A fenced code block that no line can close: the fence is longer than any run of backticks."""
    fence = "`" * max([3] + [len(run) + 1 for line in lines for run in re.findall("`+", line)])
    return [fence] + lines + [fence]


def failure_logs(work, step):
    """(log, label) of the files that tell why a step failed."""
    if step == "compare":
        return [(os.path.join(work, "output", f"compare_{m}.txt"), label) for m, label in MODES]
    if step in ("data", "mc"):
        # the output check runs only once stage1 has succeeded (an older check file is removed first)
        check = os.path.join(work, "output", f"check_{step}.txt")
        return [(check if os.path.exists(check) else os.path.join(work, "logs", f"stage1_{step}.log"), "")]
    return [(os.path.join(work, "logs", f"{step}.log"), "")]


def failures(work, who):
    """For each failed step of a work directory: its error lines, in a code block."""
    out = []
    for step in STEPS:
        if (step_status(work, step) or "").split()[:1] != ["FAIL"]:
            continue
        shown, quoted, last, last_name, empty = [], [], [], None, []
        for log, label in failure_logs(work, step):
            text = read(log)
            if text is None:
                continue
            name = f"`{os.path.join(os.path.basename(os.path.normpath(work)), os.path.relpath(log, work))}`"
            found = error_lines(text)
            if found:
                found = found[:MAX_ERROR_LINES - len(quoted)]   # at most MAX_ERROR_LINES for a step
                if found:
                    shown.append(name)
                    quoted += [f"{label}: {line}" if label else line for line in found]
                continue
            tail = last_lines(text)
            if not tail:
                empty.append(name)
            elif not last:
                last, last_name = tail, name
        if quoted:
            out += [f"**{step}**, {who}: error lines of {' and '.join(shown)}", ""] + code_block(quoted) + [""]
        elif last:
            out += [f"**{step}**, {who}: no error line found, last lines of {last_name}", ""] + code_block(last) + [""]
        elif empty:
            out += [f"**{step}**, {who}: {' and '.join(empty)} {'is' if len(empty) == 1 else 'are'} empty", ""]
        else:
            out += [f"**{step}**, {who}: no log found", ""]
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
    failed = failures(args.work, THIS) + (failures(base_work, args.base_label) if base_work else [])
    if failed:
        out += ["#### Failed steps", ""] + failed
        if os.environ.get("GITHUB_ACTIONS"):
            out += ["Full logs: job log and the logs-and-outputs artifact, on the run page.", ""]
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
