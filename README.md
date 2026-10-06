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

The key4hep stack is taken from `FCCAnalyses/.fccana/stack_pin` if that file exists, unless a key4hep stack is already set up in the shell. The file is not part of the repository, so on a fresh clone `setup.sh` sources the latest key4hep release. To use the 2025-05-29 release, create the pin once after cloning:

```bash
mkdir -p FCCAnalyses/.fccana
echo /cvmfs/sw.hsf.org/key4hep/releases/2025-05-29/x86_64-almalinux9-gcc14.2.0-opt/key4hep-stack/2025-05-30-4x4qya/setup.sh > FCCAnalyses/.fccana/stack_pin
```

Files written with the 2026-04-08 release (ROOT 6.38) are zstd-compressed and cannot be opened with the uproot of the 2025-05-29 release, so keep stage1 and the uproot-based scripts on the same release.

Set up the environment (sources the key4hep stack, configures the `FCCAnalyses` `PATH`/`PYTHONPATH`, and puts `src/` on `PYTHONPATH` for the shared modules such as `run_list`):

```bash
source setup.sh
```

Build the `FCCAnalyses` submodule (only needed once, or after pulling submodule updates or making edits to the core code yourself of course):

```bash
cd FCCAnalyses
fccanalysis build -j 8
cd ..
```

> To use another stack, write the path of its `setup.sh` to `FCCAnalyses/.fccana/stack_pin`. `FCCAnalyses/setup.sh` also accepts `-l/--latest`, `-n/--nightlies`, or `-b/--from-build` if you need a different stack than the pinned one (see `source FCCAnalyses/setup.sh --help`). After changing the stack (including creating the pin in a checkout that was already built), open a new shell, `source setup.sh` and rebuild with `fccanalysis build --clean-build -j 8`: a plain `fccanalysis build` keeps the CMake configuration of the previous stack.

Every new shell session, just re-run `source setup.sh` from the repo root before working with `fccanalysis` or the plotting scripts.

## Repository layout

- [`src/`](src/) — Stage1 (ntuple production from raw ALEPH data/MC via `fccanalysis`) and stage2 (event-level → jet-level conversion) processing. See [src/README.md](src/README.md) and [src/stage2/README.md](src/stage2/README.md).
- [`src/training/`](src/training/) — Jet-flavour tagger training configs (`weaver`).
- [`data/lumi/`](data/lumi/) — Per-year run list of the ALEPH data with the luminosity per run, and the scripts that build it from the ALEPH run database. See [src/README.md](src/README.md).
- [`Data_MC_plotting/`](Data_MC_plotting/) — Config-driven Data/MC comparison plotting for stage1 and inference-level ntuples. See [Data_MC_plotting/README.md](Data_MC_plotting/README.md).
- [`ROOT-Plotting/`](ROOT-Plotting/) — PyROOT-based plotting scripts. *Probably obsolete to be double checked*
- [`FCCAnalyses/`](FCCAnalyses/) — FCCAnalyses submodule (analyzers, build system, `fccanalysis` CLI).
- [`test/`](test/) — Validation datasets and scripts used to cross-check ntuple production across framework versions. *Probably obsolete to be double checked*

## Quick start

See [src/README.md](src/README.md) for how to run stage1 (on data or MC) and stage2, and [Data_MC_plotting/README.md](Data_MC_plotting/README.md) for making comparison plots from the resulting ntuples.