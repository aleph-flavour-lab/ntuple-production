#ifndef ALEPH_SVNEW_ANALYZERS_H
#define ALEPH_SVNEW_ANALYZERS_H

/*
  Seeded secondary-vertex finder on the secondary tracks, optionally without the
  daughters of tight V0 candidates. Lengths cm, momenta GeV.
*/

#include "ROOT/RVec.hxx"
#include "TVector3.h"
#include "edm4hep/TrackState.h"
#include "FCCAnalyses/VertexingUtils.h"
#include "FCCAnalyses/VertexFitterSimple.h"

#include <algorithm>
#include <numeric>
#include <vector>

#include "aleph_units.h"
#include "analyzer_trkaux.h"

namespace FCCAnalyses {
namespace AlephSVNew {

using ROOT::VecOps::RVec;
using AlephTrkAux::fitTracksCm;

constexpr double SVN_MPI = AlephMasses::kPiCh;

// Selection cuts
constexpr double SVN_CHI2 = 10.;        // maximum vertex chi2/ndf
constexpr double SVN_DIS_LO = 0.03;     // PV displacement window, low edge [cm]
constexpr double SVN_DIS_HI = 3.;       // PV displacement window, high edge [cm]
constexpr double SVN_SIGL_MAX = 0.10;   // maximum vertex sigma along the summed momentum [cm]
constexpr double SVN_SIGL_UNDEF = 999.; // sigma_L when undefined [cm], fails SVN_SIGL_MAX
static_assert(SVN_SIGL_UNDEF > SVN_SIGL_MAX, "an undefined sigma_L must fail the cut");
constexpr int SVN_MAX_TRK = 8;          // maximum tracks per candidate
constexpr double SVN_TRK_CHI2 = 5.;     // per-track chi2 contribution cap
constexpr double SVN_COS_POINT = 0.7;   // minimum cosPointing
constexpr double SVN_SEED_DR_MAX = 0.8; // seeds: maximum DeltaR between the two tracks
constexpr double SVN_2TRK_FSIG_MIN = 3.; // 2-track candidates: minimum 3D flight significance

// V0-track masking modes for findSVs
constexpr int SVN_MASK_NONE = 0;        // mask nothing
constexpr int SVN_MASK_MODE = 1;        // mask the daughters of tight V0 candidates

// Sigma along unit u of a packed covariance (xx, yx, yy, zx, zy, zz); SVN_SIGL_UNDEF unless var > 0
template <typename Cov>
inline double sigmaAlong(const Cov& c, const TVector3& u) {
  double x = u.x(), y = u.y(), z = u.z();
  double var = c[0] * x * x + c[2] * y * y + c[5] * z * z +
               2. * (c[1] * x * y + c[3] * x * z + c[4] * y * z);
  return (var > 0.) ? std::sqrt(var) : SVN_SIGL_UNDEF;
}

// DeltaR of (cos phi, sin phi, tanLambda), as in VertexSeed_best's pre-filter
inline double trackDeltaR(const edm4hep::TrackState& a, const edm4hep::TrackState& b) {
  return TVector3(std::cos(a.phi), std::sin(a.phi), a.tanLambda)
      .DeltaR(TVector3(std::cos(b.phi), std::sin(b.phi), b.tanLambda));
}

// u . T > 0: u in the hemisphere T points to (all u on one side if T is zero or undefined)
inline bool thrustSide(const TVector3& u, const TVector3& T) {
  return u.Dot(T) > 0.;
}

// 3D flight significance of v from PV, summed position covariances; 0 if undefined
inline double flightSig(const VertexingUtils::FCCAnalysesVertex& v,
                        const VertexingUtils::FCCAnalysesVertex& PV) {
  const TVector3 d(v.vertex.position[0] - PV.vertex.position[0],
                   v.vertex.position[1] - PV.vertex.position[1],
                   v.vertex.position[2] - PV.vertex.position[2]);
  const double L = d.Mag();
  if (!(L > 0.)) return 0.;
  double C[3][3];
  AlephTrkAux::sumCovPacked(v.vertex.covMatrix, PV.vertex.covMatrix, C);
  const TVector3 u = d.Unit();
  const double s2 = AlephTrkAux::quadFormCov(u, u, C);
  return s2 > 0. ? L / std::sqrt(s2) : 0.;
}

struct SVCand {
  VertexingUtils::FCCAnalysesVertex vtx;
  std::vector<int> trk;
  double chi2;   // normalised
  double mass;
};

// Vertex fit of np_tracks[idx]; true if it passes SVN_CHI2 and SVN_TRK_CHI2 (`group` = scratch)
inline bool svFitGroup(const RVec<edm4hep::TrackState>& np_tracks,
                       const std::vector<int>& idx, double solenoidBz,
                       RVec<edm4hep::TrackState>& group,
                       VertexingUtils::FCCAnalysesVertex& out) {
  group.clear();
  for (int k : idx) group.push_back(np_tracks[k]);
  auto v = fitTracksCm(group, np_tracks, solenoidBz);
  if ((int)v.updated_track_momentum_at_vertex.size() != (int)idx.size()) return false;
  double chi2 = v.vertex.chi2;
  if (!(chi2 == chi2) || chi2 >= SVN_CHI2) return false;
  if (v.reco_chi2.size() == idx.size())
    for (float rc : v.reco_chi2)
      if (rc > SVN_TRK_CHI2) return false;
  out = v;
  return true;
}

inline bool svPassWindows(const VertexingUtils::FCCAnalysesVertex& v,
                          const TVector3& pv, double& mass_out) {
  TVector3 x(v.vertex.position[0], v.vertex.position[1], v.vertex.position[2]);
  TVector3 d = x - pv;
  double dis = d.Mag();
  if (!(dis >= SVN_DIS_LO && dis <= SVN_DIS_HI)) return false;
  TVector3 psum(0., 0., 0.);
  double esum = 0.;
  for (const auto& tp : v.updated_track_momentum_at_vertex) {
    psum += tp;
    esum += std::sqrt(tp.Mag2() + SVN_MPI * SVN_MPI);
  }
  const double pmag = psum.Mag();
  if (pmag <= 0.) return false;
  if (!(d.Dot(psum) / (dis * pmag) >= SVN_COS_POINT)) return false;
  if (sigmaAlong(v.vertex.covMatrix, psum.Unit()) > SVN_SIGL_MAX) return false;
  double m2 = esum * esum - psum.Mag2();
  mass_out = (m2 > 0.) ? std::sqrt(m2) : 0.;
  return true;
}

// Two-track fits of all track pairs, once per event before V0 masking (findSVs filters them)
struct SVSeeds {
  std::vector<SVCand> seeds;  // window-passing pairs within SVN_SEED_DR_MAX, (i, j) ascending
  std::vector<char> pairok;   // nTr x nTr: the pair fits a common vertex
  int nTr = 0;
};

inline SVSeeds svSeedPass(const RVec<edm4hep::TrackState>& np_tracks,
                          const VertexingUtils::FCCAnalysesVertex& PV,
                          double solenoidBz) {
  SVSeeds out;
  const int nTr = np_tracks.size();
  out.nTr = nTr;
  if (nTr < 2) return out;
  TVector3 pv(PV.vertex.position[0], PV.vertex.position[1], PV.vertex.position[2]);

  RVec<edm4hep::TrackState> group;
  group.reserve(nTr);
  out.pairok.assign((size_t)nTr * nTr, 0);
  for (int i = 0; i < nTr - 1; ++i) {
    for (int j = i + 1; j < nTr; ++j) {
      VertexingUtils::FCCAnalysesVertex v;
      if (!svFitGroup(np_tracks, {i, j}, solenoidBz, group, v)) continue;
      out.pairok[(size_t)i * nTr + j] = out.pairok[(size_t)j * nTr + i] = 1;
      // growth may still attach tracks at any DeltaR
      if (trackDeltaR(np_tracks[i], np_tracks[j]) > SVN_SEED_DR_MAX) continue;
      double m;
      if (!svPassWindows(v, pv, m)) continue;
      out.seeds.push_back({v, {i, j}, v.vertex.chi2, m});
    }
  }
  return out;
}

// Secondary vertices of np_tracks (pass = its svSeedPass), V0 daughters masked per mask_mode,
// growth kept in thrustAxis hemispheres; invM = all-pion mass, reco_ind indexes np_tracks.
inline VertexingUtils::FCCAnalysesV0 findSVs(
    const RVec<edm4hep::TrackState>& np_tracks,
    const VertexingUtils::FCCAnalysesVertex& PV,
    const VertexingUtils::FCCAnalysesV0& v0s,
    const RVec<int>& v0_tight,
    int mask_mode,
    double solenoidBz,
    const SVSeeds& pass,
    const TVector3& thrustAxis) {

  VertexingUtils::FCCAnalysesV0 result;
  const int nTr = np_tracks.size();
  if (nTr < 2 || pass.nTr != nTr) return result;

  TVector3 pv(PV.vertex.position[0], PV.vertex.position[1], PV.vertex.position[2]);

  std::vector<bool> masked(nTr, false);
  if (mask_mode > 0) {
    for (size_t iv = 0; iv < v0s.vtx.size(); ++iv) {
      if (mask_mode == 1 && (iv >= v0_tight.size() || v0_tight[iv] != 1)) continue;
      for (int ti : v0s.vtx[iv].reco_ind)
        if (ti >= 0 && ti < nTr) masked[ti] = true;
    }
  }

  RVec<edm4hep::TrackState> group;
  group.reserve(nTr);
  std::vector<int> trial;
  trial.reserve(nTr);
  const std::vector<char>& pairok = pass.pairok;

  std::vector<TVector3> trkDir(nTr);
  std::vector<bool> trkSide(nTr);
  for (int k = 0; k < nTr; ++k) {
    const edm4hep::TrackState& t = np_tracks[k];
    trkDir[k] = TVector3(std::cos(t.phi), std::sin(t.phi), t.tanLambda);
    trkSide[k] = thrustSide(trkDir[k], thrustAxis);
  }

  // Grow a seed by the best-chi2 linked track in its thrust hemisphere, along the candidate
  auto growCand = [&](const SVCand& seed, const std::vector<bool>& blocked) {
    SVCand c = seed;
    const bool seedSide = thrustSide(AlephTrkAux::candMomentum(seed.vtx), thrustAxis);
    bool grew = true;
    while (grew && (int)c.trk.size() < SVN_MAX_TRK) {
      grew = false;
      SVCand best = c;
      const TVector3 candP = AlephTrkAux::candMomentum(c.vtx);
      for (int k = 0; k < nTr; ++k) {
        if (blocked[k]) continue;
        if (trkSide[k] != seedSide) continue;
        if (!(trkDir[k].Dot(candP) > 0.)) continue;
        if (std::find(c.trk.begin(), c.trk.end(), k) != c.trk.end()) continue;
        bool linked = false;
        for (int m0 : c.trk)
          if (pairok[(size_t)m0 * nTr + k]) { linked = true; break; }
        if (!linked) continue;
        trial = c.trk;
        trial.push_back(k);
        VertexingUtils::FCCAnalysesVertex v;
        if (!svFitGroup(np_tracks, trial, solenoidBz, group, v)) continue;
        double m;
        if (!svPassWindows(v, pv, m)) continue;
        if (!grew || v.vertex.chi2 < best.chi2) {
          best = {v, trial, v.vertex.chi2, m};
          grew = true;
        }
      }
      if (grew) c = best;
    }
    return c;
  };

  // seeds without a masked track, lowest chi2/ndf first
  std::vector<size_t> order;
  order.reserve(pass.seeds.size());
  for (size_t s = 0; s < pass.seeds.size(); ++s)
    if (!masked[pass.seeds[s].trk[0]] && !masked[pass.seeds[s].trk[1]])
      order.push_back(s);
  std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) {
    return pass.seeds[a].chi2 < pass.seeds[b].chi2;
  });

  std::vector<bool> used(nTr, false);
  std::vector<bool> blocked = masked;
  for (size_t s : order) {
    if (used[pass.seeds[s].trk[0]] || used[pass.seeds[s].trk[1]]) continue;
    const SVCand c = growCand(pass.seeds[s], blocked);
    for (int t : c.trk) { used[t] = true; blocked[t] = true; }
    // a 2-track candidate dropped here keeps its tracks claimed
    if (c.trk.size() == 2 && flightSig(c.vtx, PV) <= SVN_2TRK_FSIG_MIN) continue;
    result.vtx.push_back(c.vtx);
    result.pdgAbs.push_back(0);
    result.invM.push_back(c.mass);
  }
  return result;
}

inline RVec<int> candNtracks(const VertexingUtils::FCCAnalysesV0& svs) {
  RVec<int> out;
  for (const auto& v : svs.vtx) out.push_back((int)v.reco_ind.size());
  return out;
}

// SV position components relative to the PV (comp 0/1/2 = x/y/z), cm.
inline RVec<float> candDcomp(const VertexingUtils::FCCAnalysesV0& svs,
                             const VertexingUtils::FCCAnalysesVertex& PV,
                             int comp) {
  RVec<float> out;
  for (const auto& v : svs.vtx)
    out.push_back((float)(v.vertex.position[comp] - PV.vertex.position[comp]));
  return out;
}

inline RVec<float> candSigL(const VertexingUtils::FCCAnalysesV0& svs) {
  RVec<float> out;
  for (const auto& v : svs.vtx) {
    TVector3 psum(0., 0., 0.);
    for (const auto& tp : v.updated_track_momentum_at_vertex) psum += tp;
    out.push_back(psum.Mag() > 0.
                      ? (float)sigmaAlong(v.vertex.covMatrix, psum.Unit())
                      : -1.f);
  }
  return out;
}

// Per constituent track, flat over candidates: candTrkSV = candidate index, candTrkIdx =
// index in the secondary collection, candTrkOrigIdx = index in Tracks (-1 if unmapped).
inline RVec<int> candTrkSV(const VertexingUtils::FCCAnalysesV0& svs) {
  RVec<int> out;
  for (size_t i = 0; i < svs.vtx.size(); ++i)
    for (size_t k = 0; k < svs.vtx[i].reco_ind.size(); ++k) out.push_back((int)i);
  return out;
}

inline RVec<int> candTrkIdx(const VertexingUtils::FCCAnalysesV0& svs) {
  RVec<int> out;
  for (const auto& v : svs.vtx)
    for (int t : v.reco_ind) out.push_back(t);
  return out;
}

inline RVec<int> candTrkOrigIdx(const VertexingUtils::FCCAnalysesV0& svs,
                                const RVec<int>& sec2orig) {
  RVec<int> out;
  for (const auto& v : svs.vtx)
    for (int t : v.reco_ind)
      out.push_back((t >= 0 && t < (int)sec2orig.size()) ? sec2orig[t] : -1);
  return out;
}

}  // namespace AlephSVNew
}  // namespace FCCAnalyses

#endif
