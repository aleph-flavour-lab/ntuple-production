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

`--noDedxGate` (data and MC) accepts every linked dE/dx measurement as valid, i.e. switches off the failed-leg omega sentinel gate. Use it with converters that no longer copy the track omega into a failed leg.

`pfcand_d0` and `pfcand_z0` are the **raw** perigee impact parameters of the track linked to each jet constituent, taken from the stored track state: referenced to the coordinate origin, in cm, in the LCIO sign convention (the ALEPH→LCIO flip of `D0` and `omega` is applied; the `pfcand_*` covariance branches are read from the same flipped collection). They are *not* recomputed at the primary vertex — `pfcand_dxy` and `pfcand_dz` are the PV-referenced ones. A constituent with no track (a neutral) carries the guard value −9 in both, again as for the covariance branches; the PV-referenced `pfcand_dxy/dz/phi0/C/ct` also read −9 for a charge-0 constituent that carries a track link.

The −9 sentinel of `pfcand_d0`/`pfcand_z0` lies inside the physical range of those variables, so the safe test for "no track" is the track link itself (or exact equality with −9), not a range cut.

`pfcand_C` is now half the fitted curvature, |omega|/2 in 1/cm, signed by the charge. It was previously Bz·c/(2 pT) from the energy-flow pT in 1/mm, so both the scale (a factor 10) and the source differ — the latter matters for V0 daughters, whose energy-flow momentum is quoted at the V0 vertex.

`pfcand_dEdx_pads_type` and `pfcand_dEdx_wires_type` are no longer a validity mask: an accepted leg can carry any type and a rejected one reads −9 in value, error and type, so test the value branch rather than `type == 0`.

### The primary-vertex fit

The primary vertex is reconstructed by a standalone fitter, [`analyzer_pvnew.h`](analyzer_pvnew.h): a damped Gauss-Newton fit of the track helices to a common point, with a Gaussian beamspot constraint, a deterministic seed ladder, and iterative pruning of the tracks that are incompatible with the vertex. It is the default; `--oldPV` restores the previous chain.

Tracks enter the fit through a pre-selection window on the impact parameters, `|D0| < 0.75 cm` and `|Z0| < 5 cm` (`PVN_D0_MAX`, `PVN_Z0_MAX`), referenced to the run beamspot. Track/vertex compatibility is then judged at `chi2max = 5` (`PVN_CHI2_MAX`); lowering it claims fewer tracks as primary and so leaves more to the secondary finders. Both fits — the selection fit that prunes the track list and the final position fit — share one beamspot constraint, of Gaussian widths 200 µm in x, 100 µm in y and 2 cm along the beam (`PVN_BS_SIGMA_X/Y/Z`, declared in cm). Every tuned value is a named `constexpr` in that header and is not configurable from the command line.

Four `int` quality flags are written: `pv_converged` (the position fit converged), `pv_split_converged` (every pruning pass converged, not only the final fit), `pv_trivial` (fewer than two pre-selected tracks entered the fit, so the vertex carries no event information even when both fits converge) and `pv_good`, the single predicate `pv_converged && pv_split_converged && !pv_trivial` that downstream users should test. On a fit that did not converge the covariance is written as zeros while the position is still written — a nonsensical position is itself the diagnostic, and `pv_good`, not the values, is the contract. Consumers are guarded on `pv_good`: the PV-referenced jet-constituent variables fall back to the beamspot position, and the secondary-vertex and V0 collections are empty, for an event without a good PV.

The PV fit covariance `Vertex_refit_cov_xx`, `_yx`, `_yy`, `_zx`, `_zy`, `_zz` (lower-triangular, cm²) and the fit quality `Vertex_refit_chi2` (χ²/ndf) are written by both chains. Under `--oldPV` a zero covariance together with `Vertex_refit_chi2 == 0` means that no fit was run at all, because fewer than two tracks passed the pre-selection.

| flag | meaning |
| --- | --- |
| `--oldPV` | legacy PV chain, unchanged from before this module: `get_PrimaryTracks` + `VertexFitter_Tk`, with the origin-referenced `|D0| < 0.75 cm`, `|Z0| < 2 cm` pre-selection instead of the beamspot-referenced one. No `pv_*` flag branches; every other branch is bit-identical to the pre-module output, with the covariance and chi2 branches above added. Note that the beamspot constraint of the `get_PrimaryTracks` selection fit is passed in 10 µm units while its track parameters are read in cm, so that constraint is off by a factor 1000 and is effectively absent; the final `VertexFitter_Tk` fit is unaffected. Its beamspot widths and its track-compatibility cut are spelled as literals in `stage1.py` and are the same numbers as `PVN_BS_SIGMA_X/Y/Z` and `PVN_CHI2_MAX`. |

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





