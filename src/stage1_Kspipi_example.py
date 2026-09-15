"""Example on top of stage1: the two daughter tracks of every tight Ks candidate.

Runs the full stage1 chain and writes one flat block per event, ks_*, from the
tight-Ks collection Ks_event that stage1 defines (getKs(V0sNewTight_event);
getKs(V0sNew_event) would be the loose tier, tight included; Lambda_event
likewise): invariant mass and momentum of the candidate, and for each daughter
its original track index, charge, momentum at the fitted vertex, dE/dx and
(MC only) the PDG code of the linked generator particle. Same command line as
stage1, e.g.

    fccanalysis run stage1_Kspipi_example.py -- --tag test --MCflavour 5 --valid

Daughter legs are in momentum order: trk1 is the higher-momentum one. The dE/dx
value and error are -1 when the measurement is missing or invalid. The ks_* block
is standalone: its candidates are the tight Ks in their stored order, with no
index back into the v0n_* list.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stage1 import Analysis as Stage1, V0N_TRKS, DEDX_BRANCHES

KS = "Ks_event"               # the tight-Ks collection stage1 defines
KS_LEG_VARS = ("origIdx", "q", "p") + DEDX_BRANCHES


class Analysis(Stage1):

    def analyzers(self, df):
        df = super().analyzers(df)
        if not self.do_v0new:
            raise RuntimeError("stage1_Kspipi_example needs the V0 module (no --oldV0)")

        # the same cand* helpers that fill the v0n_* block work on Ks_event
        df = df.Define("n_ks",    f"int({KS}.vtx.size())")
        df = df.Define("ks_invM", f"{KS}.invM")
        df = df.Define("ks_p",    f"FCCAnalyses::AlephTruth::candP({KS})")
        # the two daughters, leg 0 = trk1 = higher momentum
        for i, t in enumerate(V0N_TRKS):
            df = df.Define(f"ks_{t}_origIdx", f"FCCAnalyses::AlephV0New::candDaughterOrigIdx({KS}, sec2origIdx, {i})")
            df = df.Define(f"ks_{t}_q",       f"FCCAnalyses::AlephV0New::candDaughterCharge({KS}, SecondaryTracks_looseBS, {i})")
            df = df.Define(f"ks_{t}_p",       f"FCCAnalyses::AlephV0New::candDaughterP({KS}, {i})")
        # dE/dx of each leg, joined through ks_trk{1,2}_origIdx like the v0n legs
        df = self._define_dedx(df, [f"ks_{t}" for t in V0N_TRKS])
        # nature of the daughters (MC only): PDG code of the generator particle
        # linked to the original track, 0 if none (MCParticles = the generator
        # collection stage1 reads)
        if self.do_truth:
            for t in V0N_TRKS:
                df = df.Define(f"ks_{t}_truePdg",
                               f"FCCAnalyses::AlephTruth::trackTruePdg(ks_{t}_origIdx, trackToMCs, MCParticles)")
        return df

    def output(self):
        legs = [f"ks_{t}_{v}" for t in V0N_TRKS for v in KS_LEG_VARS]
        if self.do_truth:
            legs += [f"ks_{t}_truePdg" for t in V0N_TRKS]
        return super().output() + ["n_ks", "ks_invM", "ks_p"] + legs
