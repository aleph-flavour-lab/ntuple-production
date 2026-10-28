"""Tight Ks -> pi+ pi- example. Finds the tight Ks candidates of stage1's V0
module and writes, per Ks, its momentum, vertex and flight, and per daughter
track its kinematics, its dE/dx with ALPHA's particle hypotheses and (on MC) its
truth. Run from this directory (src/), on data and on MC:

    fccanalysis run Kspipi_example.py -i <data.root> -o <out.root> -j 16 -- --doData
    fccanalysis run Kspipi_example.py -i <mc.root> -o <out.root> -j 16

The input is an EDM4hep file of a converter version that writes ALPHA's dE/dx
collections (dEdxWiresPID, dEdxPadsPID; on MC also dEdxWiresQDEDXMPID,
dEdxWiresQDEDXMMCLink, trackMCLink and the MCParticles parents). An input
without them (an older conversion) stops with an error naming them, before
any event is read.

Options after `--`: --doData (no truth branches, data beam spot), --MCflavour N
(MC only: keep events of one primary quark flavour, 1-5 = d u s c b, default
all), --year, --noRunList and --excludeRuns RUN [RUN ...] (the run selection of
stage1, from run_list.py; the last two data only), --beamspotJson PATH (data
only). An unknown option is an error.

The events kept are stage1's: on data the runs of the run list and event class
16, on MC the class-16 events with a primary quark (of flavour N with
--MCflavour, as stage1; of any flavour without it), and in both at least two
jets. For stage1 run with its default reconstruction options (in particular
without --oldPV and --oldTrackSel), event by event: the Ks candidate branches
and ks_trk{1,2}_{origIdx,q,p,nTPC,isChargedHad} equal stage1's v0n_* entries at
v0n_pdg == 310 && v0n_tight == 1 (ks_p is the magnitude of v0n_px, v0n_py,
v0n_pz), and the daughters' dE/dx _value, _nSamples, _expected_<h>, _sigma_<h>,
_pull_<h> and truth branches equal stage1's per-track trk_* branches at their
origIdx. ks_trk{1,2}_costheta and ks_trk{1,2}_phi have no stage1 counterpart.

Output branches. event_number, run_number and n_ks have one value per event,
every other branch one entry per tight Ks of the event. trk1 = the higher-
momentum daughter, trk2 = the other; h = ele, mu, pi, kaon, proton. Unit in [],
value when undefined in (). The dE/dx and truth branches are also -999 for a
daughter without a track (origIdx -1).
  Event
    event_number, run_number           event and run number
    n_ks                               number of tight Ks in the event
  Ks candidate
    ks_invM                            mass under the pi+ pi- hypothesis [GeV]
    ks_p, ks_px, ks_py, ks_pz          summed daughter momentum at the fitted vertex [GeV]
    ks_vx, ks_vy, ks_vz                fitted vertex position [cm]
    ks_dxyz, ks_Lxy                    3D / transverse flight length from the primary vertex [cm] (-1)
    ks_cosPointing                     cosine between flight direction and momentum (-2)
    ks_alpha                           Armenteros-Podolanski alpha, positive daughter first (-99)
    ks_qt                              Armenteros-Podolanski qT [GeV]; 0 for a vanishing Ks momentum
    ks_chi2                            vertex-fit chi2; one degree of freedom for two tracks
  Daughter kinematics
    ks_trk{1,2}_origIdx                index of the daughter track in Tracks (-1)
    ks_trk{1,2}_q                      charge [e] (0)
    ks_trk{1,2}_p                      momentum magnitude at the fitted vertex [GeV] (-1)
    ks_trk{1,2}_costheta, _phi         polar-angle cosine and azimuth [rad] of that momentum (-99)
    ks_trk{1,2}_nTPC                   TPC hits of the track (-1)
    ks_trk{1,2}_isChargedHad           1 = the track's particle is a charged hadron (energy-flow
                                       type 0, pion assumed); 0 = another type: electron or muon
                                       (flagged only above about 1 GeV; the electron identification
                                       uses dE/dx, so do not select on it in a dE/dx study) or a
                                       track of a V0 found by ALEPH's YV0V (type 3: about 10-40 % of
                                       the tight-Ks daughters, decreasing with momentum); -1 = no
                                       linked reconstructed particle of energy-flow type 0 to 5.
                                       Fractions and threshold measured on one 1994 data and one
                                       MC file.
  Daughter dE/dx, normalised to about 1 for a minimum-ionising particle; -999 = not available
    ks_trk{1,2}_dEdx_{wires,pads}_...
                                       TPC wires (dEdxWiresPID, ALPHA's QDEDX) / pads (dEdxPadsPID, QPADX):
      _value, _nSamples                measured dE/dx and its number of samples (-999)
      _expected_<h>, _sigma_<h>        ALPHA's expected dE/dx and resolution for hypothesis h,
                                       momentum error included (-999; also for a hypothesis ALPHA
                                       cannot evaluate, error code 5: low-momentum kaons and
                                       protons, mostly on the pads, even when _value is valid)
      _pull_<h>                        (value - expected_<h>) / sigma_<h> (-999, also when
                                       expected_<h> is -999 or sigma_<h> <= 0)
    ks_trk{1,2}_dEdx_wiresQDEDXM_{value,expected_<h>,sigma_<h>,pull_<h>}
                                       MC only: the same from ALPHA's parametrised simulation QDEDXM
                                       (dEdxWiresQDEDXMPID; -999, also in an event without the ALPHA
                                       MC banks)
  MC truth (PDG codes: +211 = pi+, but +11 = e-)
    ks_trk{1,2}_truePdg                ALPHA's best-match MC particle of the track (KBESTM); 0 = none,
                                       -999 = event without the ALPHA MC banks
    ks_trk{1,2}_trueMotherPdg          its first parent, 310 for the daughter of a true Ks; 0 = no
                                       best match or no parent, -999 as truePdg

About the dE/dx. On MC with a TPC simulation of version 216 or later (300 in
the MC checked for the converter's documentation) dEdxWiresPID already holds
the full-simulation dE/dx that ALPHA selects; QDEDXM is provided for
comparison. For MC runs up to 2000 (the rule of the converter's README; the
1994 MC files this example was tested on have runs below 2000) QDEDXM draws a
Gaussian around the expected dE/dx of the best-match particle, with the
constants of data run 16001, and its expected_<h> and sigma_<h> have no
momentum-error term; above run 2000 it returns the measured dE/dx. The dE/dx
branches have the names of stage1's per-track trk_dEdx_* branches, with the
prefix ks_trk{1,2}_ instead of trk_. stage1's V0-leg branches
v0n_trk{1,2}_dEdx_{wires,pads}_{value,error} are different: they read the
dEdxWires / dEdxPads collections of a dummy hypothesis. Their value equals
ks_trk{1,2}_dEdx_{wires,pads}_value where both are valid (measured on one 1994
data and one MC file: every tight-Ks daughter), but their error is ALPHA's
resolution for a particle of mass 1 GeV and charge 1, close to sigma_proton,
and the pads measurements below about 0.15-0.2 GeV in 1994 are lost there
(kept in dEdxPadsPID).

Reconstruction chain of stage1's default, every Define written out; the
constants are the ones stage1 uses, named from aleph_units.h (field) and
analyzer_pvnew.h (primary-vertex pre-selection and fit):
    event filter -> track selection -> beam spot -> primary vertex
    -> secondary tracks -> track -> particle-flow particle join
    -> V0 module, run only when the primary vertex is good (pv_good):
    V0sNew_event      = AlephV0New::findV0s(SecondaryTracks_looseBS, VertexObject_looseBS, Bz)   one fit, loose superset
    V0sNewTight_event = AlephV0New::tightV0s(V0sNew_event)   subset copy, no refit
    Ks_event          = AlephV0New::getKs(V0sNewTight_event)
    Lambda_event      = AlephV0New::getLambda(V0sNewTight_event)
All four are copies of the same fitted candidates; every cand* helper accepts
any of them (getKs(V0sNew_event) would be the loose Ks; Lambda_event the tight
Lambdas). The candidate branches use stage1's v0n_* helpers (candPcomp,
candVtxPos, candDxyz, candLxy, candCosPointing, candChi2, candAlpha, candQt).
Per daughter (leg 0 = trk1, leg 1 = trk2):
    candDaughterOrigIdx(coll, sec2origIdx, leg)             index into Tracks
    candDaughterCharge(coll, SecondaryTracks_looseBS, leg)  sign of omega of the secondary track state
    candDaughterP / candDaughterCosTheta / candDaughterPhi(coll, leg)
                                                            fitted momentum at the V0 vertex
    subdetHits(origIdx, Tracks.subdetectorHitNumbers_begin, Tracks.subdetectorHitNumbers_end,
               _Tracks_subdetectorHitNumbers, 2)            TPC hits
    legIsChargedHad(origIdx, rpOfTrack, ParticleIDNoLumiCal)  rpOfTrack: Tracks index ->
        RecoParticlesNoLumiCal index
    AlephVars::trackDedx / trackBestMatchTruth (analyzer_alephvars.h): arrays over
        Tracks that a leg reads at origIdx (ROOT::VecOps::Take)
    Other per-track quantities are reached the same way through origIdx, e.g. subdetHits
    with component 0 (VDET) or 1 (ITC) hits, or trackChi2Ndf(origIdx, Tracks.chi2,
    Tracks.ndf) (analyzer_trkaux.h).
"""
import os
import sys
from argparse import ArgumentParser
# fccanalysis loads this file by path, so this directory's modules are made importable here
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import run_list

DATA_CLASS_BIT = 16         # event class kept on data (as in stage1)
MC_FLAVOURS = {1: "dd", 2: "uu", 3: "ss", 4: "cc", 5: "bb"}   # --MCflavour values (as in stage1)
# per-run beam-spot positions of the data
BEAMSPOT_JSON = "/eos/experiment/fcc/ee/analyses/case-studies/aleph/utils/beamspot_position_data/beamspot.json"
# ALPHA dE/dx hypotheses h, in the order of AlephVars in analyzer_alephvars.h (as in stage1)
DEDX_HYPOS = ("ele", "mu", "pi", "kaon", "proton")
# converter collections the dE/dx and truth branches read; the second set MC only (as in stage1)
ALEPH_DEDX_COLLS = ("dEdxWiresPID", "_dEdxWiresPID_parameters", "dEdxPadsPID", "_dEdxPadsPID_parameters")
ALEPH_TRUTH_COLLS = ("dEdxWiresQDEDXMPID", "_dEdxWiresQDEDXMPID_parameters", "_dEdxWiresQDEDXMMCLink_from",
                     "_dEdxWiresQDEDXMMCLink_to", "_MCParticles_parents", "trackMCLink", "_trackMCLink_from",
                     "_trackMCLink_to")


class Analysis():

    def __init__(self, cmdline_args):
        parser = ArgumentParser(description='Ks example arguments (after --)')
        parser.add_argument('--doData', action='store_true',
                            help='Run on data: class filter and beam spot from the JSON, no truth branches.')
        parser.add_argument('--MCflavour', default=None, type=int,
                            help='MC only: keep events of this primary quark flavour (1 = dd ... 5 = bb); default all.')
        parser.add_argument('--year', default='1994',
                            help='Year of the run list (data).')
        parser.add_argument('--excludeRuns', nargs='+', action='extend', default=[], type=run_list.run_number, metavar='RUN',
                            help='Data only: drop these run numbers in addition to the run list.')
        parser.add_argument('--noRunList', action='store_true',
                            help='Data only: keep every run instead of the run list (--excludeRuns still applies).')
        parser.add_argument('--beamspotJson', default=os.environ.get("ALEPH_BEAMSPOT_JSON", BEAMSPOT_JSON),
                            help='Data only: per-run beam-spot positions.')
        self.ana_args, unknown = parser.parse_known_args(cmdline_args['remaining'])
        if unknown:
            print(f"----> ERROR: unrecognised arguments: {' '.join(unknown)}")
            sys.exit(1)
        if not self.ana_args.doData and (self.ana_args.excludeRuns or self.ana_args.noRunList):
            print("----> ERROR: --excludeRuns and --noRunList apply to data only (--doData).")
            sys.exit(1)
        if self.ana_args.doData and self.ana_args.MCflavour is not None:
            print("----> WARNING: Incompatible input arguments: --MCflavour defined with --doData, will be ignored.")
        if self.ana_args.MCflavour is not None and self.ana_args.MCflavour not in MC_FLAVOURS:
            print(f"----> ERROR: Requested unknown --MCflavour {self.ana_args.MCflavour}. Options: "
                  + ", ".join(f"{k} = {v}" for k, v in MC_FLAVOURS.items()))
            sys.exit(1)

        # input and output come from the fccanalysis command line: no process list
        if not cmdline_args.get('input') and not cmdline_args.get('input_file_list'):
            sys.exit("Kspipi_example.py: give the input file(s) with -i and the output file with -o")
        self.n_threads = 16      # -j on the command line overrides

        # aleph_units.h / analyzer_pvnew.h carry the constants named below
        self.include_paths = ["aleph_units.h", "analyzer.h", "analyzer_pvnew.h", "analyzer_truth.h", "analyzer_trkaux.h", "analyzer_v0new.h", "analyzer_alephvars.h"]

    def analyzers(self, df):

        # ---- converter collections of the dE/dx and truth branches ------------
        # checked before any event is read: an input of an older converter lacks them
        have = set(str(c) for c in df.GetColumnNames())
        absent = [c for c in ALEPH_DEDX_COLLS + (() if self.ana_args.doData else ALEPH_TRUTH_COLLS) if c not in have]
        if absent:
            print(f"----> ERROR: the input lacks the converter collections {', '.join(absent)}: it comes from an "
                  "older converter version; convert it again with a current one.")
            sys.exit(1)

        # ---- event filter (stage1's) ------------------------------------------
        if self.ana_args.doData:
            df = run_list.filter_runs(df, self.ana_args.year, self.ana_args.excludeRuns, self.ana_args.noRunList)
            df = df.Filter(f"AlephSelection::sel_class_filter({DATA_CLASS_BIT})(ClassBitset)")
        else:
            # class 16 with a primary quark: getJetPID is -1 otherwise
            df = df.Define("jetPID", "AlephSelection::getJetPID(ClassBitset, MCParticles)")
            df = df.Filter("jetPID > 0" if self.ana_args.MCflavour is None else f"jetPID == {self.ana_args.MCflavour}")
        # at least two jets of the exclusive two-jet clustering of the analysis particles, filtered copies of
        # RecoParticles as in stage1: RecoParticlesNoLumiCal = energy-flow types 0 to kPFMaxAnalysisType,
        # AnalysisParticles = those without the photons and neutral hadrons below kNeutralMinE;
        # RecoParticlesNoLumiCal also carries the per-leg particle-flow label below
        df = df.Define("RecoParticlesNoLumiCal_subset", "AlephSelection::filterParticles(RecoParticles, _RecoParticles_tracks.index, ParticleID, _ParticleID_particle.index, AlephSelection::isAnalysisType())")
        df = df.Define("RecoParticlesNoLumiCal",         "RecoParticlesNoLumiCal_subset.particles")
        df = df.Define("ParticleIDNoLumiCal",            "RecoParticlesNoLumiCal_subset.pid")
        df = df.Define("_RecoParticlesNoLumiCal_tracks", "RecoParticlesNoLumiCal_subset.tracks")
        df = df.Define("AnalysisParticles", "AlephSelection::filterParticles(RecoParticlesNoLumiCal_subset, AlephSelection::notSoftNeutral(AlephSelection::kNeutralMinE)).particles")
        df = df.Define("pjetc", "JetClusteringUtils::set_pseudoJets(ReconstructedParticle::get_px(AnalysisParticles), ReconstructedParticle::get_py(AnalysisParticles), "
                                "ReconstructedParticle::get_pz(AnalysisParticles), ReconstructedParticle::get_e(AnalysisParticles))")
        df = df.Define("_jet", "JetClustering::clustering_ee_kt(2, 2, 1, 0)(pjetc)")
        df = df.Filter("JetConstituentsUtils::count_jets(JetConstituentsUtils::build_constituents_cluster(AnalysisParticles, JetClusteringUtils::get_constituents(_jet))) > 1")
        df = df.Define("event_number", "EventHeader.eventNumber")
        df = df.Define("run_number", "EventHeader.runNumber")

        # ---- track selection --------------------------------------------------
        # baseline: finite track state, positive-definite covariance, chi2/ndf <= 10, TPC hits and |z0|; .tracks,
        # .trackStates and .origIdx (index into Tracks) share one order
        df = df.Define("tracks_selected_baseline_result", "AlephSelection::select_tracks_baseline(Tracks, _Tracks_trackStates, _Tracks_subdetectorHitNumbers, AlephSelection::kTrackMinTPCHits, AlephSelection::kTrackMaxAbsZ0)")
        df = df.Define("trackstates_selected_baseline", "tracks_selected_baseline_result.trackStates")
        df = df.Define("selBaselineOrigIdx", "tracks_selected_baseline_result.origIdx")

        # ---- beam spot ----------------------------------------------------------
        # MC: at the origin by construction. Data: per-run position from the
        # JSON, in 10 um units, converted to cm below.
        if self.ana_args.doData:
            df = df.Define("BeamspotVec", f'AlephSelection::get_beamspot(run_number[0], true, "{self.ana_args.beamspotJson}")')
            df = df.Define("Beamspot_x", "BeamspotVec.X()")
            df = df.Define("Beamspot_y", "BeamspotVec.Y()")
            df = df.Define("Beamspot_z", "BeamspotVec.Z()")
        else:
            df = df.Define("Beamspot_x", "0.0")
            df = df.Define("Beamspot_y", "0.0")
            df = df.Define("Beamspot_z", "0.0")
        df = df.Define("Beamspot_x_cm", "Beamspot_x*1e-3")
        df = df.Define("Beamspot_y_cm", "Beamspot_y*1e-3")
        df = df.Define("Beamspot_z_cm", "Beamspot_z*1e-3")

        # ---- primary vertex ---------------------------------------------------
        # candidates: upper bounds on the impact parameters w.r.t. the beam spot
        df = df.Define("tracks_selected_for_vertexfit_result",
                       "AlephSelection::select_tracks_impactparameters_bs(tracks_selected_baseline_result, FCCAnalyses::AlephPVNew::PVN_D0_MAX, FCCAnalyses::AlephPVNew::PVN_Z0_MAX, Beamspot_x_cm, Beamspot_y_cm, Beamspot_z_cm)")
        df = df.Define("trackstates_selected_for_vertexfit", "tracks_selected_for_vertexfit_result.trackStates")
        # the vertex fitter wants d0 and omega with the opposite sign
        df = df.Define("trackstates_selected_for_vertexfit_flipped", "AlephSelection::flipD0_copy(trackstates_selected_for_vertexfit)")
        df = df.Define("trackstates_selected_baseline_flipped", "AlephSelection::flipD0_copy(trackstates_selected_baseline)")
        # beam-spot-constrained fit that prunes the incompatible tracks;
        # pv_good (goodPV): converged, every kept track compatible, supported by at least two tracks
        df = df.Define("PVSelNew", "FCCAnalyses::AlephPVNew::select_primary_tracks(trackstates_selected_for_vertexfit_flipped, "
                                   "FCCAnalyses::AlephPVNew::beamSpot(Beamspot_x_cm, Beamspot_y_cm, Beamspot_z_cm))")
        df = df.Define("pv_good", "int(FCCAnalyses::AlephPVNew::goodPV(PVSelNew))")
        df = df.Define("RecoedPrimaryTracks_looseBS",
                       "FCCAnalyses::AlephPVNew::primaryTracksFromSel(trackstates_selected_for_vertexfit_flipped, PVSelNew, Beamspot_x_cm, Beamspot_y_cm, Beamspot_z_cm)")
        df = df.Define("VertexObject_looseBS", "FCCAnalyses::AlephPVNew::toFCCVertex(PVSelNew)")

        # ---- secondary tracks -------------------------------------------------
        # every baseline track not used by the primary vertex; sec2origIdx maps
        # them back to Tracks (truth-free track-state matching, data too)
        df = df.Define("SecondaryTracks_looseBS", "VertexFitterSimple::get_NonPrimaryTracks(trackstates_selected_baseline_flipped, RecoedPrimaryTracks_looseBS)")
        df = df.Define("sec2origIdx", "FCCAnalyses::AlephTruth::secondaryToOriginalTrack(SecondaryTracks_looseBS, trackstates_selected_baseline_flipped, selBaselineOrigIdx)")

        # ---- particle-flow join ------------------------------------------------
        # Tracks index -> index of its particle in RecoParticlesNoLumiCal (-1 = none)
        df = df.Define("rpOfTrack", "FCCAnalyses::AlephTrkAux::rpIndexByTrack(RecoParticlesNoLumiCal, _RecoParticlesNoLumiCal_tracks, Tracks.size())")

        # ---- V0 module --------------------------------------------------------
        df = df.Define("V0sNew_event",      "pv_good ? FCCAnalyses::AlephV0New::findV0s(SecondaryTracks_looseBS, VertexObject_looseBS, FCCAnalyses::AlephUnits::kBz) "
                                            ": FCCAnalyses::AlephV0New::V0Collection{}")
        df = df.Define("V0sNewTight_event", "FCCAnalyses::AlephV0New::tightV0s(V0sNew_event)")
        df = df.Define("Ks_event",          "FCCAnalyses::AlephV0New::getKs(V0sNewTight_event)")
        df = df.Define("Lambda_event",      "FCCAnalyses::AlephV0New::getLambda(V0sNewTight_event)")

        # ---- the tight Ks: mass under the pi pi hypothesis and momentum -------
        df = df.Define("n_ks",    "int(Ks_event.vtx.size())")
        df = df.Define("ks_invM", "Ks_event.invM")
        df = df.Define("ks_p",    "FCCAnalyses::AlephTruth::candP(Ks_event)")
        # summed daughter momentum at the fitted vertex [GeV]
        df = df.Define("ks_px",   "FCCAnalyses::AlephTruth::candPcomp(Ks_event, 0)")
        df = df.Define("ks_py",   "FCCAnalyses::AlephTruth::candPcomp(Ks_event, 1)")
        df = df.Define("ks_pz",   "FCCAnalyses::AlephTruth::candPcomp(Ks_event, 2)")
        # fitted vertex position [cm], flight from the primary vertex and vertex-fit chi2 (as stage1's v0n_*)
        df = df.Define("ks_vx",   "FCCAnalyses::AlephTruth::candVtxPos(Ks_event, 0)")
        df = df.Define("ks_vy",   "FCCAnalyses::AlephTruth::candVtxPos(Ks_event, 1)")
        df = df.Define("ks_vz",   "FCCAnalyses::AlephTruth::candVtxPos(Ks_event, 2)")
        df = df.Define("ks_dxyz", "FCCAnalyses::AlephTruth::candDxyz(Ks_event, VertexObject_looseBS)")
        df = df.Define("ks_Lxy",  "FCCAnalyses::AlephV0New::candLxy(Ks_event, VertexObject_looseBS)")
        df = df.Define("ks_cosPointing", "FCCAnalyses::AlephTruth::candCosPointing(Ks_event, VertexObject_looseBS)")
        df = df.Define("ks_chi2", "FCCAnalyses::AlephTruth::candChi2(Ks_event)")
        # Armenteros-Podolanski alpha (positive daughter first) and qT [GeV]
        df = df.Define("ks_alpha", "FCCAnalyses::AlephV0New::candAlpha(Ks_event, SecondaryTracks_looseBS)")
        df = df.Define("ks_qt",    "FCCAnalyses::AlephV0New::candQt(Ks_event)")

        # ---- the two daughters of every tight Ks ------------------------------
        # trk1 = leg 0 = the higher-momentum daughter, trk2 = leg 1
        df = df.Define("ks_trk1_origIdx", "FCCAnalyses::AlephV0New::candDaughterOrigIdx(Ks_event, sec2origIdx, 0)")
        df = df.Define("ks_trk2_origIdx", "FCCAnalyses::AlephV0New::candDaughterOrigIdx(Ks_event, sec2origIdx, 1)")
        df = df.Define("ks_trk1_q",       "FCCAnalyses::AlephV0New::candDaughterCharge(Ks_event, SecondaryTracks_looseBS, 0)")
        df = df.Define("ks_trk2_q",       "FCCAnalyses::AlephV0New::candDaughterCharge(Ks_event, SecondaryTracks_looseBS, 1)")
        df = df.Define("ks_trk1_p",       "FCCAnalyses::AlephV0New::candDaughterP(Ks_event, 0)")
        df = df.Define("ks_trk2_p",       "FCCAnalyses::AlephV0New::candDaughterP(Ks_event, 1)")
        # polar-angle cosine and azimuth of the same fitted momentum
        df = df.Define("ks_trk1_costheta", "FCCAnalyses::AlephV0New::candDaughterCosTheta(Ks_event, 0)")
        df = df.Define("ks_trk2_costheta", "FCCAnalyses::AlephV0New::candDaughterCosTheta(Ks_event, 1)")
        df = df.Define("ks_trk1_phi",      "FCCAnalyses::AlephV0New::candDaughterPhi(Ks_event, 0)")
        df = df.Define("ks_trk2_phi",      "FCCAnalyses::AlephV0New::candDaughterPhi(Ks_event, 1)")
        # TPC hits of the daughter track (-1 = no hit-number block)
        df = df.Define("ks_trk1_nTPC", "FCCAnalyses::AlephTrkAux::subdetHits(ks_trk1_origIdx, Tracks.subdetectorHitNumbers_begin, Tracks.subdetectorHitNumbers_end, _Tracks_subdetectorHitNumbers, 2)")
        df = df.Define("ks_trk2_nTPC", "FCCAnalyses::AlephTrkAux::subdetHits(ks_trk2_origIdx, Tracks.subdetectorHitNumbers_begin, Tracks.subdetectorHitNumbers_end, _Tracks_subdetectorHitNumbers, 2)")
        # particle-flow label: 1 = charged hadron, 0 = another type, -1 = no linked particle of type 0 to 5
        df = df.Define("ks_trk1_isChargedHad", "FCCAnalyses::AlephTrkAux::legIsChargedHad(ks_trk1_origIdx, rpOfTrack, ParticleIDNoLumiCal)")
        df = df.Define("ks_trk2_isChargedHad", "FCCAnalyses::AlephTrkAux::legIsChargedHad(ks_trk2_origIdx, rpOfTrack, ParticleIDNoLumiCal)")

        # ---- ALPHA dE/dx and best-match truth (converter collections) --------
        # AlephVars::trackDedx / trackBestMatchTruth give arrays over the whole Tracks collection,
        # as stage1's trk_dEdx_* / trk_truePdg / trk_trueMotherPdg; each leg reads its entry at origIdx.
        # -999 (AlephVars::kNotAvailable) = not available, under the conditions listed in the docstring
        AV = "FCCAnalyses::AlephVars"
        self.aleph_branches = []
        # (branch infix, as in stage1's trk_dEdx_<infix>_*; collection; parameter count; empty collection allowed)
        dets = [("wires", "dEdxWiresPID", "kNParDedx", "false"),
                ("pads", "dEdxPadsPID", "kNParDedx", "false")]
        if not self.ana_args.doData:
            dets.append(("wiresQDEDXM", "dEdxWiresQDEDXMPID", "kNParQDEDXM", "true"))
        for name, coll, npar, empty_ok in dets:
            df = df.Define(f"trkDedx_{name}", f'{AV}::trackDedx({coll}, _{coll}_parameters, Tracks.size(), '
                                              f'{AV}::{npar}, "{coll}", {empty_ok})')
        # MC only: ALPHA's best-match MC particle of each track (KBESTM) and its first parent
        if not self.ana_args.doData:
            df = df.Define("trkBestMatch", f"{AV}::trackBestMatchTruth(Tracks.size(), dEdxWiresQDEDXMPID.empty(), "
                                           "_dEdxWiresQDEDXMMCLink_from, _dEdxWiresQDEDXMMCLink_to, MCParticles, "
                                           "_MCParticles_parents, _trackMCLink_from, _trackMCLink_to, "
                                           "trackMCLink.weight)")
        for t in ("trk1", "trk2"):
            # the leg's Tracks index in the index type of Take; -1 becomes out of range (-999)
            df = df.Define(f"_ks_{t}_take", f"ROOT::VecOps::RVec<std::size_t>(ks_{t}_origIdx.begin(), ks_{t}_origIdx.end())")
            for name, _, _, _ in dets:
                members = [("value", "value")] + ([("nSamples", "nSamples")] if name != "wiresQDEDXM" else []) + \
                          [(f"{q}_{h}", f"{q}[{k}]") for k, h in enumerate(DEDX_HYPOS) for q in ("expected", "sigma")]
                for b, m in members:
                    self.aleph_branches.append(f"ks_{t}_dEdx_{name}_{b}")
                    df = df.Define(self.aleph_branches[-1], f"ROOT::VecOps::Take(trkDedx_{name}.{m}, _ks_{t}_take, {AV}::kNotAvailable)")
                for k, h in enumerate(DEDX_HYPOS):
                    self.aleph_branches.append(f"ks_{t}_dEdx_{name}_pull_{h}")
                    df = df.Define(self.aleph_branches[-1], f"ROOT::VecOps::Take(trkDedx_{name}.pull[{k}], _ks_{t}_take, {AV}::kNotAvailable)")
            if not self.ana_args.doData:
                for b, m in (("truePdg", "pdg"), ("trueMotherPdg", "motherPdg")):
                    self.aleph_branches.append(f"ks_{t}_{b}")
                    df = df.Define(self.aleph_branches[-1], f"ROOT::VecOps::Take(trkBestMatch.{m}, _ks_{t}_take, {AV}::kNotAvailableInt)")
        return df

    def output(self):
        branches = [
            "event_number", "run_number",
            "n_ks", "ks_invM", "ks_p",
            "ks_px", "ks_py", "ks_pz", "ks_vx", "ks_vy", "ks_vz", "ks_dxyz", "ks_Lxy", "ks_cosPointing", "ks_chi2", "ks_alpha", "ks_qt",
            "ks_trk1_origIdx", "ks_trk1_q", "ks_trk1_p", "ks_trk1_costheta", "ks_trk1_phi", "ks_trk1_nTPC", "ks_trk1_isChargedHad",
            "ks_trk2_origIdx", "ks_trk2_q", "ks_trk2_p", "ks_trk2_costheta", "ks_trk2_phi", "ks_trk2_nTPC", "ks_trk2_isChargedHad",
        ]
        return branches + self.aleph_branches
