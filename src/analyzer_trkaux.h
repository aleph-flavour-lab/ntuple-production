#ifndef ALEPH_TRKAUX_H
#define ALEPH_TRKAUX_H

/*
  Track-level auxiliaries and vertex-fit glue shared by the V0, phi->KK and D*
  finders and their stage1 wiring. Lengths cm, momenta GeV; the per-track
  quantities are aligned with a selected trackstate collection through its
  original-Tracks index map.
*/

#include <algorithm>
#include <cmath>
#include <vector>

#include <ROOT/RVec.hxx>
#include "TVector3.h"

#include "edm4hep/TrackState.h"
#include "edm4hep/ParticleIDData.h"
#include "FCCAnalyses/VertexingUtils.h"
#include "FCCAnalyses/VertexFitterSimple.h"

#include "aleph_units.h"

namespace FCCAnalyses {
namespace AlephTrkAux {

using ROOT::VecOps::RVec;

// Momentum at the perigee; pT = kPtPerTeslaCm*Bz/|omega| [GeV].
inline TVector3 perigeeMomentum(const edm4hep::TrackState& t, double Bz) {
  double om = std::abs(t.omega);
  if (om <= 0.) return TVector3(0., 0., 0.);
  double pt = AlephUnits::kPtPerTeslaCm * Bz / om;
  return TVector3(pt * std::cos(t.phi), pt * std::sin(t.phi), pt * t.tanLambda);
}

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

// Membership table indexed by original-Tracks index.
inline std::vector<char> memberMask(const RVec<int>& set_orig) {
  int mx = -1;
  for (int o : set_orig) mx = std::max(mx, o);
  std::vector<char> m((size_t)(mx + 1), 0);
  for (int o : set_orig)
    if (o >= 0) m[o] = 1;
  return m;
}

inline bool inMask(const std::vector<char>& m, int o) {
  return o >= 0 && (size_t)o < m.size() && m[o] != 0;
}

// 3D significance of d = x1 - x2 under the summed covariances; -1 if singular.
template <typename CovA, typename CovB>
inline float vertexDistSig(const TVector3& d, const CovA& ca, const CovB& cb) {
  double C[3][3];
  sumCovPacked(ca, cb, C);
  double det = C[0][0] * (C[1][1] * C[2][2] - C[1][2] * C[2][1])
             - C[0][1] * (C[1][0] * C[2][2] - C[1][2] * C[2][0])
             + C[0][2] * (C[1][0] * C[2][1] - C[1][1] * C[2][0]);
  if (!(std::abs(det) > 0.) || !std::isfinite(det)) return -1.;
  double inv[3][3];
  inv[0][0] = (C[1][1] * C[2][2] - C[1][2] * C[2][1]) / det;
  inv[0][1] = (C[0][2] * C[2][1] - C[0][1] * C[2][2]) / det;
  inv[0][2] = (C[0][1] * C[1][2] - C[0][2] * C[1][1]) / det;
  inv[1][0] = (C[1][2] * C[2][0] - C[1][0] * C[2][2]) / det;
  inv[1][1] = (C[0][0] * C[2][2] - C[0][2] * C[2][0]) / det;
  inv[1][2] = (C[0][2] * C[1][0] - C[0][0] * C[1][2]) / det;
  inv[2][0] = (C[1][0] * C[2][1] - C[1][1] * C[2][0]) / det;
  inv[2][1] = (C[0][1] * C[2][0] - C[0][0] * C[2][1]) / det;
  inv[2][2] = (C[0][0] * C[1][1] - C[0][1] * C[1][0]) / det;
  double dv[3] = {d.X(), d.Y(), d.Z()};
  double s2 = 0.;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) s2 += dv[i] * inv[i][j] * dv[j];
  return (s2 > 0. && std::isfinite(s2)) ? std::sqrt(s2) : -1.;
}

// Summed daughter momentum of a fitted vertex [GeV]: the candidate momentum
// every mass, pointing, flight and jet-relative quantity is built from.
inline TVector3 candMomentum(const VertexingUtils::FCCAnalysesVertex& v) {
  TVector3 p(0., 0., 0.);
  for (const auto& tp : v.updated_track_momentum_at_vertex) p += tp;
  return p;
}

// Quadratic form a^T C b of a summed position covariance (sumCovPacked output);
// a == b is the variance projected on that direction.
inline double quadFormCov(const TVector3& a, const TVector3& b,
                          const double C[3][3]) {
  const double av[3] = {a.X(), a.Y(), a.Z()}, bv[3] = {b.X(), b.Y(), b.Z()};
  double s = 0.;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) s += av[i] * C[i][j] * bv[j];
  return s;
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

// 1 for each entry whose original-Tracks index appears in set_orig.
inline RVec<int> flagInSet(const RVec<int>& orig_idx, const RVec<int>& set_orig) {
  RVec<int> out;
  const std::vector<char> in = memberMask(set_orig);
  for (int o : orig_idx) out.push_back(inMask(in, o) ? 1 : 0);
  return out;
}

// Original-Tracks indices of the daughters of the V0 candidates kept by `keep`.
inline RVec<int> claimedOrigIdx(const RVec<int>& d1, const RVec<int>& d2,
                                const RVec<int>& keep) {
  RVec<int> out;
  for (size_t i = 0; i < d1.size() && i < d2.size(); ++i) {
    if (i < keep.size() && !keep[i]) continue;
    if (d1[i] >= 0) out.push_back(d1[i]);
    if (d2[i] >= 0) out.push_back(d2[i]);
  }
  return out;
}

inline RVec<float> trackChi2Ndf(const RVec<int>& orig_idx,
                                const RVec<float>& chi2,
                                const RVec<int>& ndf) {
  RVec<float> out;
  for (int o : orig_idx) {
    float v = -1.f;
    if (o >= 0 && o < (int)chi2.size() && o < (int)ndf.size() && ndf[o] != 0)
      v = chi2[o] / float(ndf[o]);
    out.push_back(v);
  }
  return out;
}

// Per-track membership bits, in the ORIGINAL Tracks frame.
enum TrkMemberBit : int {
  kTrkPV        = 1 << 0,  // in the fitted primary-vertex track set
  kTrkV0        = 1 << 1,  // daughter of any stored V0 candidate (loose tier)
  kTrkV0Tight   = 1 << 2,  // daughter of a tight V0 candidate
  kTrkPhi       = 1 << 3,  // leg of any stored phi->KK candidate
  kTrkPhiWp     = 1 << 4,  // leg of a phi candidate passing the wp flag
  kTrkD0        = 1 << 5,  // leg of a reconstructed D0 -> K pi candidate
  kTrkDstar     = 1 << 6,  // leg of any stored D* candidate, slow pion included
  kTrkDstarTight= 1 << 7,  // leg of a D* candidate passing the tight flag
  kTrkSV        = 1 << 8,  // constituent track of a secondary vertex
  kTrkBaseline  = 1 << 9   // baseline-selected track
};

struct TrackTags {
  RVec<int> member;  // OR of TrkMemberBit, indexed by original track index
  RVec<int> nCand;   // stored V0 + phi + D0 + D* candidates using the track
};

namespace detail {

inline void tagBit(RVec<int>& m, const RVec<int>& idx, int bit) {
  for (int o : idx)
    if (o >= 0 && o < (int)m.size()) m[o] |= bit;
}

inline void tagBitIf(RVec<int>& m, const RVec<int>& idx, const RVec<int>& flag,
                     int bit) {
  for (size_t k = 0; k < idx.size() && k < flag.size(); ++k)
    if (flag[k] && idx[k] >= 0 && idx[k] < (int)m.size()) m[idx[k]] |= bit;
}

inline void countLegs(RVec<int>& n, const RVec<int>& idx) {
  for (int o : idx)
    if (o >= 0 && o < (int)n.size()) ++n[o];
}

}  // namespace detail

// Constituent track indices of every secondary vertex.
inline RVec<int> svTrackIdx(const RVec<VertexingUtils::FCCAnalysesVertex>& svs) {
  RVec<int> out;
  for (const auto& v : svs)
    for (int i : v.reco_ind) out.push_back(i);
  return out;
}

// Index lists are in the original Tracks frame, except sv_trk_idx, which
// sv2orig maps there.
// A finder that did not run passes empty lists.
inline TrackTags trackTags(size_t nTracks,
                           const RVec<int>& baseline_orig,
                           const RVec<int>& prim_orig,
                           const RVec<int>& sv_trk_idx,
                           const RVec<int>& sv2orig,
                           const RVec<int>& v0_d1, const RVec<int>& v0_d2,
                           const RVec<int>& v0_tight,
                           const RVec<int>& phi_t1, const RVec<int>& phi_t2,
                           const RVec<int>& phi_wp,
                           const RVec<int>& d0_k, const RVec<int>& d0_pi,
                           const RVec<int>& ds_k, const RVec<int>& ds_pi,
                           const RVec<int>& ds_pis, const RVec<int>& ds_tight) {
  TrackTags out;
  out.member = RVec<int>(nTracks, 0);
  out.nCand = RVec<int>(nTracks, 0);

  detail::tagBit(out.member, baseline_orig, kTrkBaseline);
  detail::tagBit(out.member, prim_orig, kTrkPV);
  for (int s : sv_trk_idx)
    if (s >= 0 && s < (int)sv2orig.size()) {
      const int o = sv2orig[s];
      if (o >= 0 && o < (int)nTracks) out.member[o] |= kTrkSV;
    }

  for (const RVec<int>* leg : {&v0_d1, &v0_d2}) {
    detail::tagBit(out.member, *leg, kTrkV0);
    detail::tagBitIf(out.member, *leg, v0_tight, kTrkV0Tight);
    detail::countLegs(out.nCand, *leg);
  }
  for (const RVec<int>* leg : {&phi_t1, &phi_t2}) {
    detail::tagBit(out.member, *leg, kTrkPhi);
    detail::tagBitIf(out.member, *leg, phi_wp, kTrkPhiWp);
    detail::countLegs(out.nCand, *leg);
  }
  for (const RVec<int>* leg : {&d0_k, &d0_pi}) {
    detail::tagBit(out.member, *leg, kTrkD0);
    detail::countLegs(out.nCand, *leg);
  }
  for (const RVec<int>* leg : {&ds_k, &ds_pi, &ds_pis}) {
    detail::tagBit(out.member, *leg, kTrkDstar);
    detail::tagBitIf(out.member, *leg, ds_tight, kTrkDstarTight);
    detail::countLegs(out.nCand, *leg);
  }
  return out;
}

// Per-daughter storage block, one instance per candidate role.
struct TrkBlock {
  RVec<int>   origIdx;   // index into the original Tracks collection
  RVec<int>   q;         // physical charge, +sign(omega)
  RVec<float> p, costheta;  // at the fitted vertex, or at the perigee [GeV]
  RVec<float> d0, z0;    // impact parameters w.r.t. the beamspot [cm]
  RVec<float> sigd0;     // sqrt(cov[0]) [cm]
  RVec<int>   nvdet, nitc;
  RVec<float> chi2ndf;   // track fit chi2/ndf
  RVec<int>   isprim;    // 1 = in the fitted primary set
  RVec<int>   pool;      // 0 = primary, 1 = secondary, 2 = neither, -1 = unstaged
};

inline void reserveTrk(TrkBlock& b, size_t n) {
  b.origIdx.reserve(n); b.q.reserve(n); b.p.reserve(n); b.costheta.reserve(n);
  b.d0.reserve(n); b.z0.reserve(n); b.sigd0.reserve(n); b.nvdet.reserve(n);
  b.nitc.reserve(n); b.chi2ndf.reserve(n); b.isprim.reserve(n); b.pool.reserve(n);
}

// Per-entry auxiliaries of one trackstate collection; a missing entry reads -1.
struct TrkAux {
  const RVec<int>*   orig_idx;
  const RVec<int>*   nvdet;
  const RVec<int>*   nitc;
  const RVec<float>* chi2ndf;
  const RVec<int>*   isprim;
  const RVec<int>*   pool;
  double bsx = 0., bsy = 0., bsz = 0.;  // beamspot [cm], reference of d0/z0

  int ival(const RVec<int>* src, int k) const {
    return (src && k >= 0 && k < (int)src->size()) ? (*src)[k] : -1;
  }
  float fval(const RVec<float>* src, int k) const {
    return (src && k >= 0 && k < (int)src->size()) ? (*src)[k] : -1.f;
  }
};

inline void pushTrk(TrkBlock& b, const RVec<edm4hep::TrackState>& tracks, int k,
                    const TVector3& p, const TrkAux& aux) {
  b.origIdx.push_back(aux.ival(aux.orig_idx, k));
  b.q.push_back((tracks[k].omega > 0) ? 1 : -1);
  b.p.push_back(p.Mag());
  b.costheta.push_back(p.Mag() > 0. ? p.Z() / p.Mag() : -99.);
  // perigee re-referenced to the beamspot; flipped states: perigee at +D0*n, n = (-sin phi, cos phi), dphi/ds = -omega
  const double cphi = std::cos(tracks[k].phi), sphi = std::sin(tracks[k].phi);
  const double s = aux.bsx * cphi + aux.bsy * sphi;
  b.d0.push_back(tracks[k].D0 + aux.bsx * sphi - aux.bsy * cphi - 0.5 * tracks[k].omega * s * s);
  b.z0.push_back(tracks[k].Z0 - aux.bsz + tracks[k].tanLambda * s);
  b.sigd0.push_back(tracks[k].covMatrix[0] > 0. ? std::sqrt(tracks[k].covMatrix[0]) : -1.);
  b.nvdet.push_back(aux.ival(aux.nvdet, k));
  b.nitc.push_back(aux.ival(aux.nitc, k));
  b.chi2ndf.push_back(aux.fval(aux.chi2ndf, k));
  b.isprim.push_back(aux.ival(aux.isprim, k));
  b.pool.push_back(aux.ival(aux.pool, k));
}

}  // namespace AlephTrkAux

// ---------------------------------------------------------------------------
// Legacy V0 finder configuration: the one source of the windows get_V0s_ALEPH
// hands to the compiled finder. Masses in GeV, displacements in cm. The loose
// tier is the wide-open booking tier; its gamma upper mass is negative, so a
// conversion is never booked.
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
