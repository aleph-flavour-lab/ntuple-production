# Aleph

Version of the FCCAnalyses code that supports command line arguments, to be able to process both data and MC with the same script. Nightlies version of the key4hep stack is required for this. 

## Setup

Folow the steps described in the general [README](../README.md) at the top level to setup the code, then just `cd src`. 
 
<!-- ```bash
git clone https://github.com/Apranikstar/Aleph.git
cd Aleph
git submodule update --init --recursive
cd FCCAnalyses
fccanalysis build -j 8
cd ../src
source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh #or compile and source the FCCAnalyses module.
```

Need nightlies because updated FCCAnalyses version after this commit is needed: https://github.com/HEP-FCC/FCCAnalyses/pull/474 -->

## Stage1: Produce ntuples

Note: Change the data fraction based on your needs.

### Run on MC:
```bash
fccanalysis run stage1.py -- --tag <version_tag>  --MCflavour <flavour_index>
```

Output files will be in: 
`/eos/experiment/fcc/ee/analyses/case-studies/aleph/processedMC/<year>/<mc_type>/stage1/<version_tag>/<flavour_name>.root`. 

Fraction of events to process can be set via `--fraction <val>`, default is to process all events. 

`--year <year>` and `--MCtype <type>` are also supported command line arguments, currently we only have `1994` and `zqq` here. 

### Run on data:
```bash
fccanalysis run stage1.py -- --tag <version_tag> --doData 
```

Output files will be in: `/eos/experiment/fcc/ee/analyses/case-studies/aleph/processedData/<year>/stage1/<version_tag>/`

`--year` and `--fraction` is also supported as an argument here. 

Only the data runs selected in `data/lumi/run_list_<year>.csv` are kept, before any selection. The list holds every run of the year in the ALEPH run database with its detector-quality flags, bookkeeping counters and stored SICAL luminosity; a run is selected when VDET, ITC and TPC are flagged P, no subsystem is flagged D, it is not in our own veto, and, when it has a stored SICAL luminosity, its bookkeeping counters cover the whole run, i.e. the hadronic event count of the files (at least 50 events) does not exceed the database Z count by more than 10 % (column `reason`). The counter check does not apply to the runs without a stored luminosity: they are valued from the lepton-pair rate of the files themselves, which the database counters do not enter. Our veto (27 runs among the flag-selected ones; column `reason` = `vetoed_*`): unusable primary vertex (VDET off, 27776 with sigma(d0) doubled and ITC hit loss), tracking-quality outliers (mean barrel TPC hit count or the fraction of tracks with sigma(d0) < 100 um more than 4 sigma below the neighbouring runs, or ITC hits per track 6-13 % below the other runs of the same fill on the raw files), the off-plateau solenoid field of 28125 (1.5176 T from the magnet-current readback and the track curvature), and three runs without dE/dx calibration banks (26420, 27750, 29667: every track carries the failed-measurement code on both the pad and the wire leg, so the analysis would have no dE/dx there). The list is built by `data/lumi/build_run_list.py` from the database and the per-run event-class counts of the converted files (`data/lumi/count_event_classes.py`). For 1994: 2056 of 2691 runs, 47.93 pb^-1.

`--excludeRuns RUN [RUN ...]` drops further runs, `--noRunList` keeps every run (only `--excludeRuns` applies); both are data-only and rejected in MC mode. The filter is an event filter, so the FCCAnalyses bookkeeping value `eventsProcessed` keeps counting the raw input. The kept runs and their luminosity are printed when the analysis graph is built (`--excludeRuns` is subtracted there, but not in `run_list.luminosity_pb`, so a production made with it needs its luminosity taken from that printout). The luminosity of the kept runs comes from the same file through `run_list.luminosity_pb(year)` (used by `Data_MC_plotting/plotting_config_stage1.py`): per run the stored SICAL luminosity, multiplied by the fill's lepton-pair factor in the fills whose SICAL bookkeeping is off (1994: fills 2122-2170, factors 1.10-1.18), and the large-angle lepton-pair rate (event class 15, calibrated on the PERF runs of July-September) for the runs without a stored luminosity; `run_list.py` also offers the stored-SICAL-only policies. The list is read where the analysis graph is built, so it has to be reachable there (`$ALEPH_RUN_LIST_<year>` overrides the path); a year without a list stops with an error unless `--noRunList` is given.

`--noDedxGate` (data and MC) accepts every linked dE/dx measurement as valid, i.e. switches off the failed-leg omega sentinel gate. Use it with converters that no longer copy the track omega into a failed leg.

The baseline track selection shared by the primary vertex fit, the V0 finder and the secondary vertex finders requires `chi2/ndf <= 10`, a finite positive-definite perigee covariance (5x5 block, Cholesky test), at least `kTrackMinTPCHits` (4) TPC hits and `|z0| <= kTrackMaxAbsZ0` (50 cm, perigee to the origin); `--oldTrackSel` (data and MC) drops the TPC-hit and `|z0|` requirements.

`pfcand_d0` and `pfcand_z0` are the **raw** perigee impact parameters of the track linked to each jet constituent, taken from the stored track state: referenced to the coordinate origin, in cm, in the LCIO sign convention (the ALEPH→LCIO flip of `D0` and `omega` is applied; the `pfcand_*` covariance branches are read from the same flipped collection). They are *not* recomputed at the primary vertex — `pfcand_dxy` and `pfcand_dz` are the PV-referenced ones. A constituent with no track (a neutral) carries the guard value −9 in both, again as for the covariance branches; the PV-referenced `pfcand_dxy/dz/phi0/C/ct` also read −9 for a charge-0 constituent that carries a track link.

The −9 sentinel of `pfcand_d0`/`pfcand_z0` lies inside the physical range of those variables, so the safe test for "no track" is the track link itself (or exact equality with −9), not a range cut.

`pfcand_C` is now half the fitted curvature, |omega|/2 in 1/cm, signed by the charge. It was previously Bz·c/(2 pT) from the energy-flow pT in 1/mm, so both the scale (a factor 10) and the source differ — the latter matters for V0 daughters, whose energy-flow momentum is quoted at the V0 vertex.

`pfcand_dEdx_pads_type` and `pfcand_dEdx_wires_type` are no longer a validity mask: an accepted leg can carry any type and a rejected one reads −9 in value, error and type, so test the value branch rather than `type == 0`.

### Run on batch:
```
fccanalysis submit stage1.py -- --tag VXX-XX --MCflavour X --batch --chunks X
```

### STAGE 2:

Don't touch stage2.py!
Open up stage2_all.py and change the desired input and output directories. 
Set the number of cpus.
Now you can decide if you want to divide each flavor into multiple files, then you can change the argument `n_final_files`.
Run it with nightlies.





