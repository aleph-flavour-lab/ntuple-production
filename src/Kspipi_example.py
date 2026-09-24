"""Tight Ks -> pi+ pi- from the V0 module: the two daughter tracks, their
momentum, dE/dx and (on MC) their nature. A standalone fccanalysis script with
the structure of stage1 (class Analysis: __init__ / analyzers / output) that
carries only the reconstruction the Ks block needs; nothing is imported from
stage1. Input and output are given on the command line (-i / -o, both
required); run from this directory, the headers are resolved relative to it:

    fccanalysis run Kspipi_example.py -i <edm4hep.root> -o <out.root> [-j 16]
    fccanalysis run Kspipi_example.py -i <edm4hep.root> -o <out.root> -- --doData

Options after `--`: --doData (no truth branches, data beam spot), --MCflavour N
(MC only: keep events of one primary quark flavour, 1-5 = d u s c b, default
all), --beamspotJson PATH (data only).

Reconstruction chain of stage1's default, every Define written out; the
constants are the ones stage1 uses, named from aleph_units.h (field) and
analyzer_pvnew.h (primary-vertex pre-selection and fit):
    event filter -> track selection -> beam spot -> primary vertex
    -> secondary tracks -> dE/dx join (wires) -> [MC] track -> generator-particle links
    -> V0 module, run only when the primary vertex is good (pv_good):
    V0sNew_event      = AlephV0New::findV0s(SecondaryTracks_looseBS, VertexObject_looseBS, Bz)   one fit, loose superset
    V0sNewTight_event = AlephV0New::tightV0s(V0sNew_event)   subset copy, no refit
    Ks_event          = AlephV0New::getKs(V0sNewTight_event)
    Lambda_event      = AlephV0New::getLambda(V0sNewTight_event)
All four are copies of the same fitted candidates; every cand* helper accepts
any of them (getKs(V0sNew_event) would be the loose Ks). Per daughter, leg 0 =
trk1 = the higher-momentum daughter, leg 1 = trk2:
    candDaughterOrigIdx(coll, sec2origIdx, leg)             index into Tracks
    candDaughterCharge(coll, SecondaryTracks_looseBS, leg)  physical charge
    candDaughterP(coll, leg)                                 |p| at the fitted vertex [GeV]
    trackQuantityByIndex(origIdx, dEdxWires.dQdx.value, dedxJoin_wires)
        wires dE/dx of the track through the per-event join dedxJoin_wires
        (Tracks index -> valid measurement index), -1 = no valid measurement;
        pads: the same with dEdxPads and a dedxJoin_pads built alike
    trackTruePdg(origIdx, trackToMCs, MCParticles)          MC only: PDG code of the
        generator particle linked to the track (211 pion, 321 kaon, 2212 proton,
        11 electron, 13 muon; PDG sign convention: +211 = pi+ but +11 = e-;
        0 = no link). A track links to every MC particle that left at least one
        hit on it, in the order the hits are met along the track; the converter
        leaves every link weight at 1.0, so the first link is taken. About 9% of
        linked tracks have several links, and for ~1% the first is not the
        particle that best matches the track momentum (mostly pi/K -> mu decays
        in flight, where the first link is usually the parent).

Output: event_number, run_number and a flat ks_* block, one entry per tight Ks:
n_ks, ks_invM, ks_p, ks_trk{1,2}_{origIdx,q,p,dEdx_wires_value,dEdx_wires_error}
and on MC ks_trk{1,2}_truePdg.
"""
import os
import sys
from argparse import ArgumentParser

DATA_CLASS_BIT = 16         # event class kept on data (as in stage1)
# per-run beam-spot positions of the data
BEAMSPOT_JSON = "/eos/experiment/fcc/ee/analyses/case-studies/aleph/utils/beamspot_position_data/beamspot.json"


class Analysis():

    def __init__(self, cmdline_args):
        parser = ArgumentParser(description='Ks example arguments (after --)')
        parser.add_argument('--doData', action='store_true',
                            help='Run on data: class filter and beam spot from the JSON, no truth branches.')
        parser.add_argument('--MCflavour', default=None, type=int,
                            help='MC only: keep events of this primary quark flavour (1 = dd ... 5 = bb); default all.')
        parser.add_argument('--beamspotJson', default=os.environ.get("ALEPH_BEAMSPOT_JSON", BEAMSPOT_JSON),
                            help='Data only: per-run beam-spot positions.')
        self.ana_args, _ = parser.parse_known_args(cmdline_args['remaining'])

        # input and output come from the fccanalysis command line: no process list
        if not cmdline_args.get('input') and not cmdline_args.get('input_file_list'):
            sys.exit("Kspipi_example.py: give the input file(s) with -i and the output file with -o")
        self.n_threads = 16      # -j on the command line overrides

        # aleph_units.h / analyzer_pvnew.h carry the constants named below
        self.include_paths = ["aleph_units.h", "analyzer.h", "analyzer_pvnew.h", "analyzer_truth.h", "analyzer_trkaux.h", "analyzer_v0new.h"]

    def analyzers(self, df):

        # ---- event filter -----------------------------------------------------
        if self.ana_args.doData:
            df = df.Filter(f"AlephSelection::sel_class_filter({DATA_CLASS_BIT})(ClassBitset)")
        elif self.ana_args.MCflavour is not None:
            df = df.Define("jetPID", "AlephSelection::getJetPID(ClassBitset, MCParticles)")
            df = df.Filter(f"jetPID == {self.ana_args.MCflavour}")
        df = df.Define("event_number", "EventHeader.eventNumber")
        df = df.Define("run_number", "EventHeader.runNumber")

        # ---- track selection --------------------------------------------------
        # baseline: positive-definite covariance, chi2 < 10, TPC hits and |z0|; .tracks,
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
        # pv_good = converged, fully pruned and supported by at least two tracks
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

        # ---- dE/dx join (wires) -----------------------------------------------
        # Tracks index -> index of its valid wires measurement (-1 = none), built
        # once per event; the last argument switches the validity gate on
        # (a failed leg carries the track omega as value)
        df = df.Define("dedxJoin_wires",
                       "FCCAnalyses::AlephV0New::dedxIndexByTrack(dEdxWires.dQdx.value, dEdxWires.dQdx.error, _dEdxWires_track.index, Tracks, _Tracks_trackStates, true)")

        # ---- MC: track -> generator-particle links ----------------------------
        if not self.ana_args.doData:
            df = df.Define("trackToMCs", "FCCAnalyses::AlephTruth::buildTrackToMCs(Tracks.size(), _trackMCLink_from, _trackMCLink_to)")

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

        # ---- the two daughters of every tight Ks ------------------------------
        # trk1 = leg 0 = the higher-momentum daughter, trk2 = leg 1
        df = df.Define("ks_trk1_origIdx", "FCCAnalyses::AlephV0New::candDaughterOrigIdx(Ks_event, sec2origIdx, 0)")
        df = df.Define("ks_trk2_origIdx", "FCCAnalyses::AlephV0New::candDaughterOrigIdx(Ks_event, sec2origIdx, 1)")
        df = df.Define("ks_trk1_q",       "FCCAnalyses::AlephV0New::candDaughterCharge(Ks_event, SecondaryTracks_looseBS, 0)")
        df = df.Define("ks_trk2_q",       "FCCAnalyses::AlephV0New::candDaughterCharge(Ks_event, SecondaryTracks_looseBS, 1)")
        df = df.Define("ks_trk1_p",       "FCCAnalyses::AlephV0New::candDaughterP(Ks_event, 0)")
        df = df.Define("ks_trk2_p",       "FCCAnalyses::AlephV0New::candDaughterP(Ks_event, 1)")
        # wires dE/dx of the daughter track, by its Tracks index through the join
        df = df.Define("ks_trk1_dEdx_wires_value", "FCCAnalyses::AlephV0New::trackQuantityByIndex(ks_trk1_origIdx, dEdxWires.dQdx.value, dedxJoin_wires)")
        df = df.Define("ks_trk1_dEdx_wires_error", "FCCAnalyses::AlephV0New::trackQuantityByIndex(ks_trk1_origIdx, dEdxWires.dQdx.error, dedxJoin_wires)")
        df = df.Define("ks_trk2_dEdx_wires_value", "FCCAnalyses::AlephV0New::trackQuantityByIndex(ks_trk2_origIdx, dEdxWires.dQdx.value, dedxJoin_wires)")
        df = df.Define("ks_trk2_dEdx_wires_error", "FCCAnalyses::AlephV0New::trackQuantityByIndex(ks_trk2_origIdx, dEdxWires.dQdx.error, dedxJoin_wires)")
        # nature of the daughter (MC only): the generator particle behind the track
        if not self.ana_args.doData:
            df = df.Define("ks_trk1_truePdg", "FCCAnalyses::AlephTruth::trackTruePdg(ks_trk1_origIdx, trackToMCs, MCParticles)")
            df = df.Define("ks_trk2_truePdg", "FCCAnalyses::AlephTruth::trackTruePdg(ks_trk2_origIdx, trackToMCs, MCParticles)")
        return df

    def output(self):
        branches = [
            "event_number", "run_number",
            "n_ks", "ks_invM", "ks_p",
            "ks_trk1_origIdx", "ks_trk1_q", "ks_trk1_p", "ks_trk1_dEdx_wires_value", "ks_trk1_dEdx_wires_error",
            "ks_trk2_origIdx", "ks_trk2_q", "ks_trk2_p", "ks_trk2_dEdx_wires_value", "ks_trk2_dEdx_wires_error",
        ]
        if not self.ana_args.doData:
            branches += ["ks_trk1_truePdg", "ks_trk2_truePdg"]
        return branches
