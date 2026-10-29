#ifndef ALEPH_RECO_CONFIG_H
#define ALEPH_RECO_CONFIG_H
// Primary-vertex and V0 pair-fit constants shared by stage1.py, the standalone
// fitter (analyzer_pvnew.h), the track helpers (analyzer_trkaux.h) and the
// standalone examples. Dependency-free;
// lengths in cm unless the name says otherwise (solenoid field and track
// conventions: aleph_units.h).
namespace FCCAnalyses {
namespace AlephReco {
// impact-parameter preselection of the primary-vertex candidate tracks
constexpr double kPVTrackD0Max = 0.75;   // |d0| [cm]
constexpr double kPVTrackZ0Max = 2.0;    // |z0| [cm]

// beam-spot constraint of the primary-vertex fit: loose widths from data
constexpr double kBeamSigmaX_um = 200.;
constexpr double kBeamSigmaY_um = 100.;
constexpr double kBeamSigmaZ_cm = 2.;
// VertexFitter_Tk scales these widths from um to mm by 1e-3, but the track
// states here are in cm, so its argument unit is 10 um. get_PrimaryTracks
// (--oldPV) scales them by the same 1e-3 but also scales the track lengths by
// 1e-3, so its constraint is 1000x too loose, i.e. effectively absent.
constexpr double kBeamSigmaXFit = kBeamSigmaX_um / 10.;
constexpr double kBeamSigmaYFit = kBeamSigmaY_um / 10.;
constexpr double kBeamSigmaZFit = kBeamSigmaZ_cm * 1e4 / 10.;

constexpr double kPVChi2Max = 5.;        // a track stays in the primary-vertex fit below this chi2

// Reach of a track's inner hits, used by the V0 pair fit: a track seen in a
// layer cannot come from a vertex outside it. Layer radii [cm] from the
// detector axis, CERN-PPE/94-170: VDET layers (Sec. 3.1), ITC outer radius
// (Sec. 3.2).
constexpr double kVdetInnerR = 6.5;
constexpr double kVdetOuterR = 11.3;
constexpr double kItcOuterR = 26.;
// Tolerances [cm] added to those radii: 99 % quantile of the fitted minus the
// true transverse vertex radius of simulated K0S and Lambda decaying at 4-8,
// 8-14 and 22-32 cm; they also cover the flat VDET wafers.
constexpr double kHitReachTolVdetInner = 1.1;
constexpr double kHitReachTolVdetOuter = 1.8;
constexpr double kHitReachTolItc = 2.9;
}  // namespace AlephReco
}  // namespace FCCAnalyses
#endif
