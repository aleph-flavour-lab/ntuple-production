# Aleph ntuple production for flavour tagging 

Analysis framework for reviving and reprocessing archival **ALEPH** (LEP) data and MC using the [key4hep](https://key4hep.github.io/key4hep-doc/) stack and [FCCAnalyses](https://github.com/Apranikstar/FCCAnalyses). The repo covers the ntuple production step following the staged analyses approach in `FCCAnalyses`, as well as scripts for training the tagger. *Note: Code structure to be tiedied up in future.*

Since custom changes in FCCAnalyses core code are needed, a specific FCCAnalyses submodule, pinned to a fork ([Apranikstar/FCCAnalyses](https://github.com/Apranikstar/FCCAnalyses)), is required and included in this repo. 


## Setup

Clone the repo together with the `FCCAnalyses` submodule:

```bash
git clone https://github.com/aleph-flavour-lab/ntuple-production.git
cd ntuple-production
git submodule update --init --recursive
```

Set up the environment (sources the pinned key4hep stack recorded in `FCCAnalyses/.fccana/stack_pin`, configures the `FCCAnalyses` `PATH`/`PYTHONPATH`, and puts `src/` on `PYTHONPATH` for the shared modules such as `run_list`):

```bash
source setup.sh
```

Build the `FCCAnalyses` submodule (only needed once, or after pulling submodule updates or making edits to the core code yourself of course):

```bash
cd FCCAnalyses
fccanalysis build -j 8
cd ..
```

> To update the pinned stack, edit `FCCAnalyses/.fccana/stack_pin`. `FCCAnalyses/setup.sh` also accepts `-l/--latest`, `-n/--nightlies`, or `-b/--from-build` if you need a different stack than the pinned one (see `source FCCAnalyses/setup.sh --help`).

Every new shell session, just re-run `source setup.sh` from the repo root before working with `fccanalysis` or the plotting scripts.

## Repository layout

- [`src/`](src/) — Stage1 (ntuple production from raw ALEPH data/MC via `fccanalysis`) and stage2 (event-level → jet-level conversion) processing. See [src/README.md](src/README.md) and [src/stage2/README.md](src/stage2/README.md).
- [`src/training/`](src/training/) — Jet-flavour tagger training configs (`weaver`).
- [`data/lumi/`](data/lumi/) — Per-year run list of the ALEPH data with the luminosity per run, and the scripts that build it from the ALEPH run database. See [src/README.md](src/README.md).
- [`Data_MC_plotting/`](Data_MC_plotting/) — Config-driven Data/MC comparison plotting for stage1 and inference-level ntuples. See [Data_MC_plotting/README.md](Data_MC_plotting/README.md).
- [`ROOT-Plotting/`](ROOT-Plotting/) — PyROOT-based plotting scripts. *Probably obsolete to be double checked*
- [`FCCAnalyses/`](FCCAnalyses/) — FCCAnalyses submodule (analyzers, build system, `fccanalysis` CLI).
- [`test/`](test/) — Validation datasets and scripts used to cross-check ntuple production across framework versions. *Probably obsolete to be double checked*
- [`.github/`](.github/) — CI workflows and the scripts they run. See [Continuous integration](#continuous-integration).

## Quick start

See [src/README.md](src/README.md) for how to run stage1 (on data or MC) and stage2, and [Data_MC_plotting/README.md](Data_MC_plotting/README.md) for making comparison plots from the resulting ntuples.

## Continuous integration

Every pull request to `main` runs stage1 (`--doData` and `--MCflavour 5`) with GitHub Actions on two small synthetic files (data-like and MC-like toy Z → qq̄ events with the collections of the converted ALEPH files, generated in the job), for the pull request and for `main`, and compares the two outputs event by event. Pushes to `main` and manual runs (Actions tab) run stage1 without the comparison.

- The check **fails** if the pull request breaks the run (FCCAnalyses build, input generation, stage1, or an empty stage1 output) or if the comparison cannot be made. Output differences are **informational**; a failure on the `main` side only gives a warning.
- The job summary on the run page, also posted as a comment on the pull request, shows the status of each step (with the error lines of a failed one) and the branches added, removed and changed. Logs and outputs: artifact `stage1-ci-logs-and-outputs` of the run, kept 7 days.
- The synthetic events exercise the code, not the physics: identical outputs do not prove identical outputs on real data, and code that the synthetic events do not reach is not tested.

The CI has its own key4hep stack, pinned in `.github/ci/key4hep_stack` (not `FCCAnalyses/.fccana/stack_pin`). To run the same steps by hand (AlmaLinux 9 with `/cvmfs`, a new shell without a key4hep setup; `-h` lists the steps and options):

```bash
.github/ci/run_stage1_ci.sh
```

If `FCCAnalyses` was already built with another stack, the build step stops: use a separate clone.
