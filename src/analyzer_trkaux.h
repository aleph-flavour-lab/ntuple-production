#ifndef ALEPH_TRKAUX_H
#define ALEPH_TRKAUX_H

/*
  Track-level auxiliaries and vertex-fit glue shared by the V0 module and its
  stage1 wiring. Lengths cm, momenta GeV; the per-track quantities are aligned
  with a selected trackstate collection through its original-Tracks index map.
*/

#include <algorithm>
#include <cmath>

#include <ROOT/RVec.hxx>
#include "TVector3.h"

#include "edm4hep/TrackState.h"
#include "edm4hep/ParticleIDData.h"
#include "FCCAnalyses/VertexingUtils.h"
#include "FCCAnalyses/VertexFitterSimple.h"

namespace FCCAnalyses {
namespace AlephTrkAux {

using ROOT::VecOps::RVec;

// The one vertex-fit entry point of every finder: fit `group` (alltracks passed
// so reco_ind is filled), then undo the cm-as-mm homothety of the momenta,
// which come out 10x too small, ONCE at the source.
inline VertexingUtils::FCCAnalysesVertex fitTracksCm(
    const RVec<edm4hep::TrackState>& group,
    const RVec<edm4hep::TrackState>& alltracks, double Bz) {
  auto v = VertexFitterSimple::VertexFitter_Tk(
      0, group, alltracks, false, 0., 0., 0., 0., 0., 0., Bz, false);
  for (auto& tp : v.updated_track_momentum_at_vertex) tp *= 10.;
  return v;
}

// Undefined Armenteros-Podolanski variables: the value the candidate-level
// alpha/qt carry when the pair is unusable, far below any physical alpha.
constexpr double kApUndef = -99.;

// Smallest track multiplicity of a usable primary vertex: below it the PV is
// the default object at the origin, so distances measured from it are meaningless.
constexpr int kPVMinTracks = 2;
// Undefined value of a float branch (missing input, guarded quantity).
constexpr float kUndef = -1.f;

// Armenteros-Podolanski variables of a track pair: qt and
// alpha = (pL+ - pL-)/(pL+ + pL-). q1sign = physical charge of p1; for a
// same-charge pair the labels are conventional, so pass +1 to order by p1.
// alpha_null is returned when the longitudinal momenta cancel.
inline void apVars(const TVector3& p1, const TVector3& p2, double q1sign,
                   double& alpha, double& qt, double alpha_null = 0.) {
  TVector3 p = p1 + p2;
  double pmag = p.Mag();
  qt = p1.Cross(p.Unit()).Mag();
  double la = p1.Dot(p) / pmag, lb = p2.Dot(p) / pmag;
  double lplus = (q1sign > 0) ? la : lb, lminus = (q1sign > 0) ? lb : la;
  alpha = (lplus + lminus != 0.) ? (lplus - lminus) / (lplus + lminus)
                                 : alpha_null;
}

// Sum of two packed lower-triangular position covariances (xx, yx, yy, zx, zy,
// zz) as a full 3x3 matrix.
template <typename CovA, typename CovB>
inline void sumCovPacked(const CovA& ca, const CovB& cb, double C[3][3]) {
  C[0][0] = double(ca[0]) + cb[0];
  C[0][1] = C[1][0] = double(ca[1]) + cb[1];
  C[1][1] = double(ca[2]) + cb[2];
  C[0][2] = C[2][0] = double(ca[3]) + cb[3];
  C[1][2] = C[2][1] = double(ca[4]) + cb[4];
  C[2][2] = double(ca[5]) + cb[5];
}

// ---------------------------------------------------------------------------
// Particle-flow join: original track index -> ReconstructedParticle -> PF type.
// ---------------------------------------------------------------------------

// PF type code of the ParticleID collection (index-parallel to RecoParticles);
// the same source pfcand_isChargedHad reads.
constexpr int kPFChargedHad = 0;

// ReconstructedParticle index of every original-Tracks index, -1 when no RP
// points at the track. tracks_begin/tracks_end are offsets into the RP->Track
// relation, and begin != end is the "has a track" test; a neutral's begin is
// in range but meaningless.
template <typename Idx>
inline RVec<int> rpIndexByTrack(const RVec<Idx>& tracks_begin,
                                const RVec<Idx>& tracks_end,
                                const RVec<int>& rpTrackIndex, size_t nTracks) {
  RVec<int> out(nTracks, -1);
  const size_t n = std::min(tracks_begin.size(), tracks_end.size());
  for (size_t i = 0; i < n; ++i)
    for (size_t k = tracks_begin[i]; k < (size_t)tracks_end[i]; ++k) {
      if (k >= rpTrackIndex.size()) break;
      const int t = rpTrackIndex[k];
      if (t >= 0 && (size_t)t < nTracks && out[t] < 0) out[t] = (int)i;
    }
  return out;
}

// Tri-state PF charged-hadron flag of a candidate leg: 1 = its linked
// ReconstructedParticle is a PF charged hadron, 0 = the RP is e/mu/other,
// -1 = the track has no linked RP.
inline RVec<int> legIsChargedHad(const RVec<int>& orig_idx,
                                 const RVec<int>& rp_of_track,
                                 const RVec<edm4hep::ParticleIDData>& pid) {
  RVec<int> out;
  for (int o : orig_idx) {
    int v = -1;
    if (o >= 0 && o < (int)rp_of_track.size()) {
      const int r = rp_of_track[o];
      if (r >= 0 && r < (int)pid.size())
        v = (pid[r].type == kPFChargedHad) ? 1 : 0;
    }
    out.push_back(v);
  }
  return out;
}

// Component `which` (0 = VDET, 1 = ITC, 2 = TPC) of the per-track
// subdetectorHitNumbers block, by original-Tracks index. -1 when the block is
// missing or too short. Same source the per-constituent hit counts read.
inline RVec<int> subdetHits(const RVec<int>& orig_idx,
                            const RVec<unsigned int>& begin,
                            const RVec<unsigned int>& end,
                            const RVec<int>& values, int which) {
  RVec<int> out;
  for (int o : orig_idx) {
    int v = -1;
    if (o >= 0 && o < (int)begin.size() && o < (int)end.size()) {
      unsigned int b = begin[o], e = end[o];
      if (b + which < e && b + which < values.size()) v = values[b + which];
    }
    out.push_back(v);
  }
  return out;
}

}  // namespace AlephTrkAux

// ---------------------------------------------------------------------------
// Legacy V0 finder configuration: the ONE source of the windows handed to the
// compiled finder and of the replica that recovers its pair indices. Masses in
// GeV, displacements in cm. The loose tier is the wide-open ML/booking tier;
// its gamma upper mass is negative, so a conversion is never booked.
// ---------------------------------------------------------------------------
namespace AlephLegacyV0 {

constexpr double kLooseKsMLo    = 0.1,  kLooseKsMHi    = 1.4;
constexpr double kLooseLamMLo   = 0.1,  kLooseLamMHi   = 1.4;
constexpr double kLooseGammaMLo = 0.0,  kLooseGammaMHi = -1.;
constexpr double kLooseCosKs = 0.999, kLooseCosLam = 0.999, kLooseCosGamma = 0.999;

constexpr double kTightKsMLo    = 0.453, kTightKsMHi    = 0.553;
constexpr double kTightLamMLo   = 1.06,  kTightLamMHi   = 1.16;
constexpr double kTightGammaMLo = 0.0,   kTightGammaMHi = 0.005;
constexpr double kTightCosKs = 0.999, kTightCosLam = 0.99995, kTightCosGamma = 0.99995;

// common to both tiers
constexpr double kDisMinKs = 0.1, kDisMinLam = 0.1, kDisMinGamma = 0.9;
constexpr double kChi2Cut = 10.;

}  // namespace AlephLegacyV0
}  // namespace FCCAnalyses

#endif  // ALEPH_TRKAUX_H
