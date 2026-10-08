#!/usr/bin/env bash
# stage1 CI on synthetic input: the entry point of .github/workflows/ci.yml, and
# runnable by hand from any clone on an AlmaLinux 9 machine with /cvmfs (or inside the
# ghcr.io/key4hep/key4hep-images/alma9 image), from a shell without a Key4hep setup:
#
#   .github/ci/run_stage1_ci.sh [-j THREADS] [-w WORKDIR] [-n EVENTS] [-s SEED]
#                               [-i INPUTDIR] [-r REF_WORKDIR [-R RTOL] [-A ATOL]] [STEP ...]
#
# -j: threads of the FCCAnalyses build and of each stage1 run, processes of the input
#     generation and of the comparison. -n: events per file (data and MC).
# -R, -A: relative and absolute tolerance of the comparison (default 0: exact).
# Steps (default: build inputs data mc, in this order):
#   build    build FCCAnalyses, unless a build made with the CI stack is already there
#   inputs   write the synthetic data and MC files and their beam-spot file (make_synthetic_input.py)
#   data     stage1 --doData on the synthetic data, then check the output
#   mc       stage1 --MCflavour 5 on the synthetic MC, then check the output
#   compare  compare the data and MC outputs with those in REF_WORKDIR (-r), event by
#            event; differences are reported, not failed on (see step_compare)
# data and mc run at the same time when both are requested (stage1 spends most of its
# time compiling the analysis, on one core). Every requested step runs even if an
# earlier one failed; the exit status is the number of failed steps (0 = all passed).
# Logs and outputs: WORKDIR/{logs,output} (default WORKDIR: $TMPDIR/stage1_ci_$USER).
# The synthetic input is WORKDIR/input, or INPUTDIR (-i) to reuse the input of another run.
# The Key4hep stack is read from .github/ci/key4hep_stack.

CI_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO=$(cd "$CI_DIR/../.." && pwd)
NTHREADS=4   # a GitHub runner has 4 cores
WORK="${TMPDIR:-/tmp}/stage1_ci_$(id -un)"
NEVENTS=1500
SEED=1
INPUT=
REF=
RTOL=0
ATOL=0
STAGE1_TIMEOUT=1800
export PYTHONDONTWRITEBYTECODE=1   # no __pycache__ in the checkout

usage() { sed -n '2,24p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }
while getopts "j:w:n:s:i:r:R:A:h" opt; do
  case $opt in
    j) NTHREADS=$OPTARG ;;
    w) WORK=$OPTARG ;;
    n) NEVENTS=$OPTARG ;;
    s) SEED=$OPTARG ;;
    i) INPUT=$OPTARG ;;
    r) REF=$OPTARG ;;
    R) RTOL=$OPTARG ;;
    A) ATOL=$OPTARG ;;
    h) usage 0 ;;
    *) usage 1 ;;
  esac
done
shift $((OPTIND - 1))
STEPS=("$@")
[[ ${#STEPS[@]} -eq 0 ]] && STEPS=(build inputs data mc)
for s in "${STEPS[@]}"; do
  [[ $s =~ ^(build|inputs|data|mc|compare)$ ]] || { echo "unknown step '$s'"; usage 1; }
  [[ $s == compare && -z $REF ]] && { echo "the compare step needs -r REF_WORKDIR"; usage 1; }
  [[ $s == inputs && -n $INPUT ]] && { echo "-i INPUTDIR reuses an input: no inputs step"; usage 1; }
done

mkdir -p "$WORK/logs" "$WORK/output"
WORK=$(cd "$WORK" && pwd)
INPUT=${INPUT:-$WORK/input}
STACK=$(grep -v -e '^[[:space:]]*#' -e '^[[:space:]]*$' "$CI_DIR/key4hep_stack" | head -1)
SUMMARY=()
NFAIL=0

annotate() {  # annotate LEVEL TITLE MESSAGE: an annotation in GitHub Actions, plain text elsewhere
  if [[ -n ${GITHUB_ACTIONS:-} ]]; then echo "::$1 title=$2::$3"; else echo "${1^^} [$2] $3"; fi
}

fail_early() {  # fail_early LOG TITLE MESSAGE: MESSAGE as the step's log (for job_summary.py) and an annotation
  echo "ERROR: $3" > "$1"
  annotate error "$2" "$3"
}

setup_failed() {  # setup_failed MESSAGE: the setup log, then MESSAGE (also added to the log, for job_summary.py)
  cat "$WORK/logs/setup.log"
  echo "ERROR: $1" >> "$WORK/logs/setup.log"
  annotate error setup "$1"
}

setup_env() {
  local log=$WORK/logs/setup.log
  if [[ -n ${KEY4HEP_STACK:-} && $KEY4HEP_STACK != "$STACK" ]] || [[ -n ${FCCANA_LOCAL_DIR:-} && $FCCANA_LOCAL_DIR != "$REPO/FCCAnalyses" ]]; then
    fail_early "$log" setup "this shell already has another Key4hep/FCCAnalyses setup; run from a fresh shell"
    return 1
  fi
  set --   # the setup scripts parse the positional parameters
  # shellcheck disable=SC1090
  source "$STACK" > "$log" 2>&1
  # stop here if the stack did not set itself up: FCCAnalyses/setup.sh would otherwise source
  # another stack (its pinned one, or the latest release) without saying so
  if [[ ${KEY4HEP_STACK:-} != "$STACK" ]]; then
    setup_failed "the Key4hep stack of .github/ci/key4hep_stack did not set up (KEY4HEP_STACK='${KEY4HEP_STACK:-}', expected $STACK)"
    return 1
  fi
  # shellcheck disable=SC1091
  source "$REPO/setup.sh" >> "$log" 2>&1
  if [[ ${KEY4HEP_STACK:-} != "$STACK" || ${FCCANA_LOCAL_DIR:-} != "$REPO/FCCAnalyses" ]]; then
    setup_failed "setup.sh did not set up $REPO/FCCAnalyses with the CI stack (FCCANA_LOCAL_DIR='${FCCANA_LOCAL_DIR:-}', KEY4HEP_STACK='${KEY4HEP_STACK:-}')"
    return 1
  fi
  echo "stack: $KEY4HEP_STACK"
  echo "repo:  $REPO ($(git -C "$REPO" log -1 --format='%h %s' 2>/dev/null))"
  echo "FCCAnalyses: $(git -C "$REPO/FCCAnalyses" log -1 --format='%h' 2>/dev/null)"
}

step_build() {
  local lib=$REPO/FCCAnalyses/install/lib/libFCCAnalyses.so
  local built_with
  built_with=$(cat "$REPO/FCCAnalyses/.fccana/stack_build" 2>/dev/null)
  if [[ -f $lib && $built_with == "$KEY4HEP_STACK" ]]; then
    echo "FCCAnalyses already built with this stack (CI cache, or an earlier build): not rebuilding"
    return 0
  fi
  if [[ -e $REPO/FCCAnalyses/install || -e $REPO/FCCAnalyses/build ]]; then
    fail_early "$WORK/logs/build.log" build "FCCAnalyses has a build for another stack (${built_with:-unknown}); remove FCCAnalyses/build and FCCAnalyses/install, or rebuild with $STACK"
    return 1
  fi
  echo "building FCCAnalyses with $NTHREADS threads (log: $WORK/logs/build.log)"
  (cd "$REPO/FCCAnalyses" && fccanalysis build -j "$NTHREADS") > "$WORK/logs/build.log" 2>&1
  local rc=$?
  hash -r
  if [[ $rc -ne 0 || ! -f $lib ]]; then
    tail -30 "$WORK/logs/build.log"
    annotate error build "fccanalysis build failed (exit status $rc)"
    return 1
  fi
}

step_inputs() {
  python "$CI_DIR/make_synthetic_input.py" --out "$INPUT" --nevents "$NEVENTS" --seed "$SEED" \
    --jobs "$NTHREADS" > "$WORK/logs/inputs.log" 2>&1
  local rc=$?
  cat "$WORK/logs/inputs.log"
  [[ $rc -eq 0 ]] || { annotate error inputs "synthetic input generation failed (exit status $rc)"; return 1; }
}

output_name() {  # the stage1 output file of a mode, relative to WORK/output
  # names given by src/stage1.py (data: process 1994; MC: output name of --MCflavour 5)
  case $1 in
    data) echo data/1994.root ;;
    mc)   echo mc/Zbb.root ;;
  esac
}

step_stage1() {
  local mode=$1 args
  case $mode in
    data) args=(--doData) ;;
    mc)   args=(--MCflavour 5) ;;
  esac
  local root log=$WORK/logs/stage1_$mode.log
  root=$WORK/output/$(output_name "$mode")
  # stale results of an earlier run in the same WORKDIR must not pass the check (nor be
  # taken by job_summary.py as the reason of a failure)
  rm -f "$root" "$root.failed_check" "$WORK/output/branches_$mode.txt" "$WORK/output/check_$mode.txt"
  if [[ ! -f $INPUT/beamspot.json ]]; then
    fail_early "$log" "stage1 $mode" "no synthetic input in $INPUT: run the inputs step first"
    return 1
  fi
  mkdir -p "$(dirname "$root")"
  echo "stage1 ${args[*]} with $NTHREADS threads (log: $log)"
  CI_INPUT_DIR=$INPUT CI_OUTPUT_DIR=$(dirname "$root") CI_NTHREADS=$NTHREADS \
    ALEPH_BEAMSPOT_JSON=$INPUT/beamspot.json \
    timeout "$STAGE1_TIMEOUT" fccanalysis run "$CI_DIR/stage1_ci.py" -- --tag ci "${args[@]}" > "$log" 2>&1
  local rc=$?
  [[ $rc -eq 124 ]] && echo "ERROR: stage1 stopped at the time limit of $STAGE1_TIMEOUT s" >> "$log"
  if [[ $rc -ne 0 ]]; then
    local first
    first=$(grep -m1 -E 'error:|ERROR|Error in|Traceback|terminate called' "$log")
    echo "---- stage1 $mode failed (exit status $rc); error lines of the log:"
    grep -E 'error:|ERROR|Error in|runtime_error|terminate called|Traceback' "$log" | head -20
    echo "---- last lines of the log:"
    tail -15 "$log"
    annotate error "stage1 $mode" "exit status $rc: ${first:-see the log}"
    rm -f "$root"   # a partial output (crash, time limit) must not be compared
    return 1
  fi
  python "$CI_DIR/check_stage1_output.py" "$root" --branches-out "$WORK/output/branches_$mode.txt" \
    > "$WORK/output/check_$mode.txt" 2>&1
  rc=$?
  cat "$WORK/output/check_$mode.txt"
  [[ $rc -eq 0 ]] && return 0
  annotate error "stage1 $mode" "output check failed"
  # out of the comparison (an existing output counts as a good one there), kept in the artifact
  [[ -f $root ]] && mv -f "$root" "$root.failed_check"
  return 1
}

# Compares this run's outputs (WORK) with the reference run's (REF; in the CI: main).
# POLICY (informational): differences are listed in compare_<mode>.txt/.json and in the job
# summary; the step fails only when a comparison cannot be made although both outputs exist.
# A missing output (stage1 or its check failed there) gives "not compared", reported by its own step.
step_compare() {
  local mode out ref new rc nfail=0
  for mode in data mc; do
    out=$WORK/output/compare_$mode
    rm -f "$out.txt" "$out.json"
    ref=$REF/output/$(output_name "$mode")
    new=$WORK/output/$(output_name "$mode")
    if [[ ! -f $new ]]; then
      echo "not compared: no $mode output of this run (stage1 failed or did not run)" | tee "$out.txt"
      continue
    fi
    if [[ ! -f $ref ]]; then
      echo "not compared: no $mode output of the reference (stage1 failed or did not run there)" | tee "$out.txt"
      annotate warning "compare $mode" "no $mode output of the reference run: nothing to compare with"
      continue
    fi
    echo "---- $mode: reference $ref"
    python "$CI_DIR/compare_stage1_outputs.py" "$ref" "$new" --rtol "$RTOL" --atol "$ATOL" \
      --jobs "$NTHREADS" --json "$out.json" > "$out.txt" 2>&1
    rc=$?
    cat "$out.txt"
    if [[ $rc -ge 2 ]]; then
      annotate error "compare $mode" "the comparison could not be made: $(grep -m1 'comparison not possible' "$out.txt")"
      nfail=$((nfail + 1))
    fi
  done
  [[ $nfail -eq 0 ]]
}

# the commit under test, for the summary: the plain short hash, except for GitHub's test merge
# of a pull request (GITHUB_SHA in a pull_request run), shown with the pull request's head
# commit (second parent) and main (first parent); nothing if git cannot read the repository
commit_label() {
  local sha short p1 p2 more
  read -r sha short p1 p2 more < <(git -C "$REPO" log -1 --format='%H %h %p' 2>/dev/null)
  if [[ ${GITHUB_EVENT_NAME:-} == pull_request && -n $sha && $sha == "${GITHUB_SHA:-}" && -n $p2 && -z $more ]]; then
    echo "test merge $short of $p2 into $p1"
  else
    echo "$short"
  fi
}

record() {  # record STEP RC SECONDS: the summary line, also left in WORK/logs/STEP.status
  local status=PASS
  [[ $2 -ne 0 ]] && { status=FAIL; NFAIL=$((NFAIL + 1)); }
  SUMMARY+=("$(printf '%-7s %-5s %5ss' "$1" "$status" "$3")")
  echo "${SUMMARY[-1]}" > "$WORK/logs/$1.status"   # read by job_summary.py
}

run_step() {
  local step=$1 t0 rc
  echo "==================== $step ===================="
  t0=$(date +%s)
  case $step in
    build) step_build ;;
    inputs) step_inputs ;;
    data|mc) step_stage1 "$step" ;;
    compare) step_compare ;;
  esac
  rc=$?
  record "$step" "$rc" $(( $(date +%s) - t0 ))
}

# data and mc at the same time; their output is shown one after the other
run_stage1_both() {
  local mode
  rm -f "$WORK"/logs/stage1_{data,mc}.rc
  for mode in data mc; do
    ( t0=$(date +%s)
      step_stage1 "$mode" > "$WORK/logs/stage1_$mode.console" 2>&1
      rc=$?
      echo "$rc $(( $(date +%s) - t0 ))" > "$WORK/logs/stage1_$mode.rc"
    ) &
  done
  wait
  for mode in data mc; do
    echo "==================== $mode ===================="
    cat "$WORK/logs/stage1_$mode.console"
    local rc_t
    read -r -a rc_t 2>/dev/null < "$WORK/logs/stage1_$mode.rc" || rc_t=(1 0)   # no file: the run died
    record "$mode" "${rc_t[0]}" "${rc_t[1]}"
  done
}

T_START=$(date +%s)
echo "==================== setup ===================="
rm -f "$WORK/logs/setup.status"
if ! setup_env; then
  echo "stage1 CI: environment setup FAILED, no step run"
  record setup 1 $(( $(date +%s) - T_START ))   # its log, logs/setup.log, says why
  STEPS=()   # the summary below is still written
fi
both=0   # 1: data and mc both requested, run them at the same time (2: done)
[[ " ${STEPS[*]} " == *" data "* && " ${STEPS[*]} " == *" mc "* ]] && both=1
for step in "${STEPS[@]}"; do
  if [[ $both -ne 0 && $step =~ ^(data|mc)$ ]]; then
    [[ $both -eq 1 ]] && run_stage1_both
    both=2
  else
    run_step "$step"
  fi
done

# data and MC are meant to have the same branches; a difference is reported, not failed
if [[ -f $WORK/output/branches_data.txt && -f $WORK/output/branches_mc.txt ]]; then
  if ! cmp -s "$WORK/output/branches_data.txt" "$WORK/output/branches_mc.txt"; then
    echo "WARNING: data and MC outputs have different branches:"
    diff "$WORK/output/branches_data.txt" "$WORK/output/branches_mc.txt" | grep '^[<>]' | sed 's/^</  data only:/; s/^>/  MC only:  /'
  fi
fi

echo "==================== summary ===================="
{
  label=$(commit_label)
  echo "stage1 CI on synthetic input (${label:+$label, }$NTHREADS threads, $NEVENTS events per file)"
  printf '%s\n' "${SUMMARY[@]}"
  echo "total   $(( $(date +%s) - T_START ))s, failed steps: $NFAIL"
} | tee "$WORK/summary.txt"
exit "$NFAIL"
