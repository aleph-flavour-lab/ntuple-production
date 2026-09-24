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

### Reconstruction modules: defaults and opt-outs

A flag-less `stage1.py` runs the two-tier V0 module ([`analyzer_v0new.h`](analyzer_v0new.h)) on top of the standard chain, writing the `v0n_*` branches. The module writes no generator-level information, so the output schema is the same on MC and on data.

The legacy code path remains available as an opt-out:

| flag | meaning |
| --- | --- |
| `--oldV0` | drop the two-tier V0 module: no `v0n_*` branches. The legacy `v0_*` block is unaffected. |
| `--noV0TagVars` | drop the jet-relative tagger inputs of the module (the `v0n_jetIdx`/`z`/`zL`/... block and the per-daughter `q`/`p`/`nTPC` branches). The rest of the module is unaffected; `--oldV0` drops them too. |

The selection itself is not configurable from the command line: every tuned value is a named constant in the header that uses it, declared once. The paragraphs below describe it in words.

### The two-tier V0 module

Standalone V0 (Ks/Λ) reconstruction in [`analyzer_v0new.h`](analyzer_v0new.h). It runs on the primary vertex of the standard chain (`VertexObject_looseBS`) and on the secondary tracks that vertex leaves unclaimed (`SecondaryTracks_looseBS`). It is *V0-first*: its track claims are meant to be consumed by later finders, so it optimises the correctness of each claim rather than only the candidate list.

**Candidate building.** All opposite-charge secondary track pairs are vertexed with a single consistent fit (`VertexFitter_Tk`); every downstream quantity (momenta, invariant masses under both hypotheses, Armenteros–Podolanski (AP) variables, pointing) is derived from the refitted momenta at that vertex — there is no second fit.

**Two selection tiers, evaluated per hypothesis (Ks and Λ):**

- **Tight** — the adopted physics selection: mass window; momentum-tiered pointing cut (separate Ks and Λ ladders, `ksPointThr` / `lamPointThr`); a qT veto against photon conversions (Λ only); and a resolution-scaled AP-band cut around the exact kinematic locus — the band half-width follows the measured σ_ell(p) of each species (`ksBandThr`, `lamBandThrTight`; the Λ band is floored at low p and capped at the nominal ramp edge), plus common fit-quality (χ²) and displacement requirements.
- **Loose** — the ML-training tier: same windows/χ²/displacement, but flat pointing, a wider Ks AP band, a Λ AP band equal to a fixed fraction (0.8) of the ramp half-width floored at the tight band, and a relaxed Λ qT veto. It is a strict superset of tight and is what gets stored, so *any* tighter selection can be re-derived offline from any production.

**Hypothesis arbitration.** A pair passing both hypotheses is booked as the one whose invariant mass is closer to the PDG mass of that species, each distance measured in units of that species' mass-window half-width; an exact tie is booked as a Ks.

**Exclusive claiming, tight first.** Candidates claim their tracks exclusively in quality order: all tight candidates claim before any loose one, and within a tier the best-χ² candidate claims first. A track is claimed once; later candidates using it are dropped. This preserves the tight-only output exactly regardless of the loose tier.

**Stored flags and ML inputs.** `v0n_tight` is the tier the finder booked the candidate in (booking a candidate ≠ selecting it — the loose tier is stored with the flag off). The tight tier and the two hypotheses can be taken as collections of their own with the header's `tightV0s` / `getKs` / `getLambda` (copies of the same fitted candidates; `Kspipi_example.py` shows the chain). `v0n_bandSig` and `v0n_massSig` store the AP-band and mass cut variables as signed pulls in resolution units for training; all other cut variables (cosPointing, pointSig, qT, χ², displacement, p, invM) are stored raw.

**Output branches.** One block, written by default and dropped by `--oldV0`: `n_v0n_event`, `v0n_pdg`, `v0n_invM`, `v0n_alpha`, `v0n_qt`, `v0n_chi2`, `v0n_dxyz`, `v0n_px/py/pz`, `v0n_cosPointing`, `v0n_pointSig`, `v0n_tight`, `v0n_bandSig`, `v0n_massSig`, `v0n_vx/vy/vz`, the vertex-fit covariance `v0n_cov_*`, the daughter joins `v0n_trk1_origIdx`/`v0n_trk2_origIdx` and their `v0n_trk{1,2}_dEdx_{pads,wires}_{value,error}` and `v0n_trk{1,2}_isChargedHad` — event-order candidate quantities, independent of jet assignment. The candidate momentum magnitude is not written: it is \|(`v0n_px`, `v0n_py`, `v0n_pz`)\|.

Per candidate, in the order of `v0n_pdg`:

| branch | definition | unit | undefined |
| --- | --- | --- | --- |
| `v0n_pdg` | \|PDG\| of the booked hypothesis, 310 (Ks) or 3122 (Λ or Λ̄; the sign is in `v0n_alpha` / `v0n_baryon`) | — | — |
| `v0n_invM` | invariant mass under the booked hypothesis, from the refitted momenta | GeV | — |
| `v0n_alpha` | Armenteros–Podolanski α, positive track first (see below) | — | −99 |
| `v0n_qt` | Armenteros–Podolanski qT, daughter momentum transverse to the candidate momentum | GeV | 0 for a vanishing candidate momentum |
| `v0n_chi2` | vertex-fit χ² per degree of freedom (one degree of freedom for a two-track vertex) | — | — |
| `v0n_dxyz` | 3D flight length of the candidate vertex from the primary vertex | cm | −1 (primary vertex with fewer than 2 tracks) |
| `v0n_px/py/pz` | summed daughter momentum at the fitted vertex | GeV | — |
| `v0n_cosPointing` | cosine between the flight direction from the primary vertex and the candidate momentum | — | −2 (zero flight length or momentum) |
| `v0n_pointSig` | transverse pointing significance: the candidate-to-primary-vertex vector projected on the plane transverse to the candidate momentum, in units of the summed candidate and primary position covariances (√(dᵀC⁻¹d) in that plane) | — | −1 (degenerate geometry or singular covariance) |
| `v0n_tight` | 1 = booked in the tight tier, 0 = loose tier | — | — |
| `v0n_bandSig` | signed AP-band pull of the booked hypothesis, (bandEll − 1)/σ_ell(p) | — | −999 |
| `v0n_massSig` | signed mass pull of the booked hypothesis, (invM − m_hyp)/σ_m(p) | — | −999 |
| `v0n_vx/vy/vz` | fitted vertex position | cm | — |
| `v0n_cov_{xx,yx,yy,zx,zy,zz}` | vertex-fit position covariance, packed lower triangle | cm² | — |
| `v0n_trk{1,2}_origIdx` | index of the daughter in the `Tracks` collection | — | −1 |
| `v0n_trk{1,2}_dEdx_{pads,wires}_{value,error}` | daughter dE/dx measurement after the validity gate | as `pfcand_dEdx_*` | −1 |
| `v0n_trk{1,2}_isChargedHad` | tri-state particle-flow label of the daughter (see below) | — | −1 (no linked reconstructed particle) |

The two daughter legs `v0n_trk1_*` / `v0n_trk2_*` are in momentum order, `v0n_trk1_*` the higher-momentum daughter at the fitted vertex, *not* in charge order.

`v0n_alpha` follows the physical charge — the positive track is taken first, so α = (p∥⁺ − p∥⁻)/(p∥⁺ + p∥⁻), and α > 0 means Λ rather than Λ̄.

A daughter leg whose dE/dx measurement is missing or fails the validity gate reads −1 in both `v0n_trk{1,2}_dEdx_{pads,wires}_value` and `..._error`, not the −9 of the `pfcand_dEdx_*` block. `--noDedxGate` switches the gate off for both blocks alike.

### Jet-relative V0 tagger inputs

A second block of branches, written by default and dropped by `--noV0TagVars` (and by `--oldV0`), presents each candidate as the jet sees it. They are derived from the stored candidates, the primary vertex `VertexObject_looseBS` and the jet collection — no candidate is re-fitted and no tuned value enters. Every candidate is assigned to the jet with the smallest ΔR between its momentum and the jet axis, the first jet winning a tie; a candidate with a vanishing momentum is left unassigned. It is the same rule the legacy `v0_*` block is assigned by.

Momentum fractions are normalised to the energy of the assigned jet, taken from the jet collection the clustering produced. Every float branch of this block uses −1 as the undefined value; the reliable "no jet" test is `v0n_jetIdx == -1`.

Per candidate, in the order of `v0n_pdg`:

| branch | definition | unit | undefined |
| --- | --- | --- | --- |
| `v0n_jetIdx` | index of the assigned jet | — | −1 (no jet) |
| `v0n_z` | \|p_V0\| / E_jet | — | −1 |
| `v0n_zL` | (p_V0 · ĵ) / E_jet, ĵ the jet direction | — | −1 |
| `v0n_ptRel` | \|p_V0 × ĵ\|, momentum transverse to the jet axis | GeV | −1 |
| `v0n_dRjet` | ΔR between the candidate momentum and the jet axis | — | −1 |
| `v0n_rankInJet` | momentum rank among the candidates of the same jet, 1 = leading; an exact tie is broken by candidate order | — | −1 |
| `v0n_Lxy` | transverse flight length of the candidate vertex from the primary vertex (the 3D one is `v0n_dxyz`) | cm | −1 |
| `v0n_LxySig` | `v0n_Lxy` over its uncertainty: candidate and primary position covariances summed and projected on the transverse flight direction | — | −1 |
| `v0n_LxyzSig` | the same in 3D, i.e. `v0n_dxyz` over the uncertainty projected on the 3D flight direction | — | −1 |
| `v0n_baryon` | +1 Λ, −1 Λ̄, 0 for a Ks or an undefined candidate; the sign of `v0n_alpha`, i.e. the charge of the proton leg | — | 0 |
| `v0n_nShared` | daughters of this candidate (0, 1 or 2) whose original track is also a daughter of another stored candidate; 0 everywhere while claiming is exclusive | — | — |
| `v0n_trk{1,2}_q` | physical charge of the daughter | e | 0 (no daughter) |
| `v0n_trk{1,2}_p` | daughter momentum magnitude at the fitted vertex | GeV | −1 |
| `v0n_trk{1,2}_nTPC` | TPC hits of the daughter's original track, from the same per-track block as `pfcand_nTrackHits_TPC` | — | −1 |

The `v0n_trk{1,2}_*` legs of this group follow the same momentum order. `v0n_trk{1,2}_q` is the only branch giving the sign of a leg, so under `--noV0TagVars` a leg is identified by `v0n_trk{1,2}_origIdx` alone.

`v0n_Lxy`, `v0n_LxySig`, `v0n_LxyzSig` and `v0n_dxyz` are −1 when the primary vertex has fewer than 2 tracks: that vertex is the default one at the origin, so there is no flight to measure. The pointing branches (`v0n_cosPointing`, `v0n_pointSig`) and the legacy `v0_dxyz` keep measuring from that default vertex.

[`Kspipi_example.py`](Kspipi_example.py) is a standalone fccanalysis script with the structure of `stage1.py` (own `class Analysis`, nothing imported from stage1) that carries only the reconstruction the tight-Kₛ block needs: event filter, track selection, beam-spot-constrained primary vertex, secondary tracks and their `sec2origIdx` map, the wires dE/dx join, on MC the track → generator links, then `findV0s` → `tightV0s` → `getKs` / `getLambda`. It writes `event_number`, `run_number` and a flat `ks_*` block (`n_ks`, `ks_invM`, `ks_p`, per leg `ks_trk{1,2}_{origIdx,q,p,dEdx_wires_value,dEdx_wires_error}` via `candDaughterOrigIdx` / `candDaughterCharge` / `candDaughterP` / `trackQuantityByIndex`, and on MC `ks_trk{1,2}_truePdg` via `AlephTruth::trackTruePdg`, the PDG code of the generator particle linked to the daughter track; a track links to every generator particle that left a hit on it, in hit order and all with weight 1.0, and the first link is taken, so for a leg that decayed in flight (about 2%) the label is usually the parent). Run it from `src/` with `fccanalysis run Kspipi_example.py -i <edm4hep.root> -o <out.root> [-- --doData]`; its `ks_*` entries equal the `v0n_*` entries of stage1 at `v0n_pdg == 310 && v0n_tight == 1`.

**Further utilities in the headers.** [`analyzer_truth.h`](analyzer_truth.h) holds the truth-free secondary-track index recovery `secondaryToOriginalTrack`, the event-order candidate accessors (`candChi2`, `candDxyz`, `candP`, `candPcomp`, `candCosPointing`, `candVtxPos`) and the track→MC link map and per-track generator PDG lookup that `Kspipi_example.py` uses. [`analyzer_trkaux.h`](analyzer_trkaux.h) holds the shared vertex-fit glue, the track→`ReconstructedParticle` join behind the per-leg particle-flow label and the per-track subdetector hit-count lookup. Both are loaded unconditionally: everything stage1 takes from them is truth-free and runs on data.

### Index map and per-leg labels

Internally, `sec2origIdx` maps the secondary track collection back to the original `Tracks` index space, which is the frame the `v0n_trk{1,2}_origIdx` branches use; the matching is truth-free, so it runs on data too. The map itself is not an output branch. Every candidate leg additionally carries `<leg>_isChargedHad`, a tri-state particle-flow label: 1 = the leg's linked reconstructed particle is a PF charged hadron, 0 = it is an electron, muon or another type, −1 = the track has no linked reconstructed particle (mostly very soft tracks, below the PF momentum reach). The label uses the same particle-flow type codes as `pfcand_isChargedHad`.

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





