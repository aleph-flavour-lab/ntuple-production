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

### Reconstruction modules: defaults and opt-outs

A flag-less `stage1.py` runs the two-tier V0 module ([`analyzer_v0new.h`](analyzer_v0new.h)), the φ→K⁺K⁻ finder ([`analyzer_phikk.h`](analyzer_phikk.h)) and the D*→D⁰π finder ([`analyzer_dstar.h`](analyzer_dstar.h)) on top of the standard chain, writing the `v0n_*`, `phikk_*` and `dstar_*` branches. None of them writes generator-level information, so the output schema is the same on MC and on data.

The legacy code paths remain available as opt-outs:

| flag | meaning |
| --- | --- |
| `--oldV0` | drop the two-tier V0 module: no `v0n_*` branches. The legacy `v0_*` block is unaffected. The φ→KK and D* finders still run, with their Kₛ/Λ track veto empty (that veto list is a V0-module product). |
| `--noV0TagVars` | drop the jet-relative tagger inputs of the module (the `v0n_jetIdx`/`z`/`zL`/... block and the per-daughter `q`/`p`/`nTPC` branches). The rest of the module is unaffected; `--oldV0` drops them too. |
| `--noPhiKK` | skip the φ→K⁺K⁻ finder: no `phikk_*` branches and no `n_phikk_event`. The φ bits of `trk_member` then stay 0. |
| `--noDstar` | skip the D*→D⁰π finder: no `dstar_*` branches and no `n_dstar_event`/`n_d0fits_event`. The D⁰/D* bits of `trk_member` then stay 0. |

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

### The φ→K⁺K⁻ finder

Standalone φ(1020)→K⁺K⁻ reconstruction in [`analyzer_phikk.h`](analyzer_phikk.h), on by default and skipped with `--noPhiKK` (`phikk_*` branches). It is an extension of the standalone V0 machinery, intended to deliver a **kinematically tagged kaon sample for dE/dx calibration** — so **no dE/dx quantity enters any selection**; the daughters' dE/dx measurements are stored, never cut on.

**Candidate building.** Track pairs are formed from the *full* baseline-selected track list — primary and secondary tracks alike, with no masking by the PV split, minus the tight-claimed V0 daughters (see below) — and vertexed with the same single consistent `VertexFitter_Tk` call as the V0 module (momenta rescaled once by the cm-as-mm factor 10). A pre-fit K⁺K⁻ mass window on the perigee momenta removes the bulk of the pair combinatorics before any fit. Same-charge pairs are reconstructed too and flagged (`phikk_same_sign`), as the data-driven combinatorial control. There is **no exclusive claiming**: a track may appear in several candidates.

**V0 veto.** Tracks already claimed as daughters of a tight Kₛ/Λ candidate are removed from the pairing pool, the same claim list the D* finder uses. It is a product of the V0 module, so under `--oldV0` the veto is simply empty and the finder still runs. Candidates still claim nothing themselves: within the surviving pool a track may appear in several φ candidates.

**Promptness is not required.** The φ has cτ ≈ 46 fm, but φ from b/c decays (B→φX, D_s→φπ) are genuinely displaced and must be kept, so no displacement window and no pointing cut are applied. |vtx − PV| and its 3D significance are stored (`phikk_dpv`, `phikk_dpvSig`) and never cut on. A *storage fiducial* `DPV_FID` bounds the stored sample: near-collinear K⁺K⁻ pairs leave the vertex unconstrained along the flight direction, so a genuine φ can be reconstructed far from the production point.

**Armenteros–Podolanski for equal masses.** With equal daughter masses the AP locus is the ellipse (α/α_max)² + (q_T/p\*)² = 1 centred at α = 0, with p\* = √(m_φ²/4 − m_K²) and α_max = p\*/(β E\*), E\* = m_φ/2. `phikk_bandEll` stores the left-hand side (1 on the exact locus). Note that for equal masses (α, q_T) together with the pair momentum determine the invariant mass exactly, so the AP band is a momentum-dependent reparametrisation of the mass window rather than independent information; it is stored, not cut on.

**Selection.** Not configurable from the command line: every value is a named `constexpr` in `analyzer_phikk.h`, namespace `AlephPhiKK`, read directly by `findPhiKK`, so `stage1.py` passes none of them; the numbers are those declarations. The stored sample is deliberately loose, so that any working point is re-derivable offline: a K⁺K⁻ mass window (`M_LO`/`M_HI`, widened for the pre-fit stage by `PRE_MARGIN`), a loose vertex χ² ceiling (`CHI2_CUT`, sanity only, ndf = 1), a storage fiducial on |vtx − PV| (`DPV_FID`), and a daughter momentum floor (`P_MIN_DEF`, applied to the *perigee* momentum while the stored `phikk_trk*_p` is the at-vertex one). No displacement, pointing, AP-band or track-quality cut is applied, and same-charge pairs are kept as the control sample.

**Working-point flags.** Two per-candidate integer flags are stored alongside the loose sample, evaluated on the stored quantities themselves (post-fit mass, at-vertex daughter momenta): `phikk_wp` = |m − m_φ| < `WP_DM` and both daughters p > `WP_PDAU` and |vtx − PV| < `WP_DPV`; `phikk_tight` = the same with |m − m_φ| < `TIGHT_DM` and both daughters σ(d0) < `TIGHT_SIGD0`. Both are **charge-blind** on purpose: the signal sample is `phikk_wp && !phikk_same_sign` and the same-charge control sample under *identical* cuts is `phikk_wp && phikk_same_sign` (likewise for `phikk_tight`), which is what makes the sideband subtraction well defined. They are labels, not cuts — no candidate is dropped by them.

**Stored per-daughter block** (`phikk_trk1_`, `phikk_trk2_`; `trk1` is the higher-momentum leg). Original-track index, charge, momentum, cosθ, d0/z0 w.r.t. the run beamspot, σ(d0) at the perigee, nVDET/nITC hits, track χ²/ndf, whether the track was in the fitted primary set, the pads/wires dE/dx value and error, and the particle-flow label `isChargedHad` — everything the offline quality and purity scans need.

### The D*→D⁰π finder

Standalone D*⁺→D⁰(K⁻π⁺)π⁺_slow reconstruction in [`analyzer_dstar.h`](analyzer_dstar.h), on by default and skipped with `--noDstar` (`dstar_*` branches). Like the φ→KK finder it is an extension of the standalone V0 machinery whose purpose is a **kinematically tagged kaon sample for dE/dx calibration** — so **no dE/dx quantity enters any selection**; the daughters' dE/dx measurements are stored, never cut on. The mass difference Δm = m(Kππ_s) − m(Kπ) is the handle: its resolution is set by the slow pion alone, so a narrow Δm window together with the D⁰ mass window isolates a sample in which the track given the kaon mass really is a kaon.

**One output collection.** Only the D* list is written (`dstar_*`, one entry per D⁰ candidate × third track); `n_dstar_event` counts it. The D⁰→Kπ candidates the D* is built from stay internal, and `dstar_d0idx` is the index of the parent entry in that internal list: entries sharing a `dstar_d0idx` share the same K/π pair, mass assignment and fitted vertex, and differ only by the slow pion. `n_d0fits_event` counts the two-track fits actually performed, i.e. the combinatorial cost of the event.

**Candidate building.** Track pairs are formed from the *full* baseline-selected track list — primary and secondary tracks alike, with no masking by the PV split — and vertexed with the same single consistent `VertexFitter_Tk` call as the V0 module (momenta rescaled once by the cm-as-mm factor 10). **Both mass assignments** of every opposite-charge pair are separate candidates: the kaon hypothesis is what defines the tag, so (a=K, b=π) and (a=π, b=K) are different objects, not a symmetry to be resolved. A pre-fit Kπ mass window on the perigee momenta removes the bulk of the pair combinatorics before any fit, and a pair is fitted once however many mass assignments and slow pions reach it. There is **no exclusive claiming**: a track may appear in several candidates. Per candidate, `dstar_nsec` records how many legs sit in the secondary pool, and each daughter carries `..._pool` (0 primary set, 1 secondary set, 2 neither) beside `..._isprim`.

**No three-track fit.** The D⁰ flies ≈0.6 mm while the slow pion comes from the D* decay point, i.e. from the PV region, so a common three-track vertex would be wrong. The D* is built from the D⁰ momenta **at the fitted D⁰ vertex** plus the slow pion's **perigee** momentum, and Δm is required below `DM_MAX`.

**Right-sign vs wrong-sign.** The slow pion of a true D*⁺ carries the charge of the pion from the D⁰ (both opposite to the kaon). Both slow-pion charges are kept and flagged (`dstar_rs` = 1 for right-sign). The wrong-sign sample is **not** signal-free and is not a drop-in background template: the flag compares the slow pion with the track *given the pion mass*, so every true D* is also reconstructed with the K↔π masses swapped, and that duplicate is wrong-sign by construction. Δm is nearly insensitive to the K/π assignment, so the duplicate peaks at the **same** Δm as the signal: in simulation studies it is concentrated inside the `DS_LOOSE_DDM` window and, before any D⁰ mass window, amounts to roughly 40 % of the true signal. The genuinely combinatorial part of the wrong-sign sample *is* flat in Δm, so it is this reflection alone that breaks the estimator: a wrong-sign template normalised in a Δm **sideband** is normalised where the reflection is absent and therefore over-subtracts inside the peak, removing a sizeable fraction of the signal and biasing the tagged-kaon dE/dx distribution the finder exists to provide. The same deficit is visible in data, as a dip of the right-sign/wrong-sign ratio confined to the D⁰-mass region where the reflection lands.

Offline remedies, in increasing order of effort: tighten the D⁰ mass window, since the swapped assignment is displaced in m(Kπ) even though Δm is not; or veto wrong-sign candidates whose swapped partner also falls in the D⁰ window — both mass assignments of a pair are separate entries sharing the same daughter `*_origIdx`, so the partner is identifiable in data as well as in simulation; or normalise the wrong-sign template to the in-peak background level rather than to a Δm sideband. Neither a plain RS − WS difference nor a sideband-normalised one is unbiased, and neither should be used for an absolute yield without one of these corrections.

**Promptness is not required.** D* from b decays give a D⁰ that does not point at the PV, so there is no displacement window and no pointing cut: `dstar_dpv` and `dstar_dpvSig` (3D significance vs the PV) are **stored, never cut on**, and only the wide storage fiducial `DPV_FID` bounds |vtx − PV|; `dstar_cosPoint` is stored too and enters the `dstar_tight` label alone, which rejects random pairs without requiring promptness. `dstar_cosThetaStar`, the cosine of the kaon direction in the D⁰ rest frame w.r.t. the D⁰ lab flight direction, is stored too: flat for a true two-body decay, peaked at |cos| = 1 for combinatorics.

**Selection.** Not configurable from the command line: every value is a named `constexpr` in `analyzer_dstar.h`, namespace `AlephDstar`, read directly by `findDstar`, so `stage1.py` passes none of them; the numbers are those declarations. What is *stored* is deliberately loose, so any working point is re-derivable offline: a Kπ mass window (`M_LO`/`M_HI`, widened for the pre-fit stage by `PRE_MARGIN`), a loose D⁰ vertex χ² ceiling (`CHI2_CUT`, sanity only), a storage fiducial on |vtx − PV| (`DPV_FID`), a Δm ceiling (`DM_MAX`), and perigee momentum floors for the K and π (`P_MIN`) and for the slow pion (`PS_MIN`). No track-quality prefilter is applied.

**V0 veto.** Tracks already claimed as daughters of a tight Kₛ/Λ candidate are removed from the pool before any pairing. The claim list is a product of the V0 module, so under `--oldV0` the veto is simply empty and the finder still runs.

**Working-point flags.** Labels, not cuts — no candidate is dropped by them, and they are evaluated on the stored post-fit quantities. `dstar_loose` = |m(Kπ) − m_D⁰| < `DS_LOOSE_DM` and |Δm − `DM_NOMINAL`| < `DS_LOOSE_DDM`, with no requirement on the slow pion or the pointing. `dstar_tight` tightens those to `DS_TIGHT_DM` and `DS_TIGHT_DDM` and adds p(K) > `TIGHT_PK`, p(π) > `TIGHT_PPI`, a D⁰ vertex χ² below `TIGHT_CHI2`, p(π_s) > `DS_TIGHT_PS`, `cosPoint` > `DS_TIGHT_COSPOINT`, and a **primary-pattern veto**: a candidate whose K and π are both in the fitted primary set while its slow pion is not is rejected, that pattern being combinatorial rather than a real D* topology.

**Stored per-daughter block** (`dstar_trkK_`, `dstar_trkPi_`, and `dstar_trkPis_` for the slow pion). Original-track index, charge, momentum, cosθ, d0/z0 (beamspot-referenced, as on the φ legs), σ(d0) at the perigee, nVDET/nITC hits, track χ²/ndf, whether the track was in the fitted primary set, the staging pool, the pads/wires dE/dx value and error, and the particle-flow label `isChargedHad`.

### Track membership and joins

`pfcand_trackIdx`, `trk_member` and `trk_nCand` are written unconditionally — with both finders off and under `--oldV0` as well. The bits and the candidate counts of a finder that did not run simply stay 0.

`pfcand_trackIdx` gives the **original `Tracks`** index of each jet constituent's own track (−1 for a constituent with no track). That is the index space every finder's `*_origIdx` branch already uses, so the `pfcand_*` block and any candidate leg can be matched directly instead of through the `ReconstructedParticle`→`Track` relation.

`trk_member` and `trk_nCand` are per-track arrays over the whole `Tracks` collection. `trk_member` is a bitmask recording every set a track belongs to: bit 0 fitted primary-vertex set, bit 1 daughter of any stored V0 candidate, bit 2 daughter of a *tight* V0, bit 3 leg of any stored φ→KK candidate, bit 4 leg of a φ candidate passing `phikk_wp`, bit 5 leg of a reconstructed D⁰→Kπ candidate, bit 6 leg of any stored D* candidate (slow pion included), bit 7 leg of a D* passing `dstar_tight`, bit 8 constituent track of a secondary vertex found by `get_SV_event_ALEPH`, bit 9 baseline-selected track. `trk_nCand` counts how many candidates (V0 + φ + D⁰ + D*) use the track, which is the multiplicity the offline 1/n de-duplication weight needs — the finders claim no tracks exclusively, so one track can serve many candidates. Both are built in a single pass over the finished candidate lists, so they add no reconstruction work.

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





