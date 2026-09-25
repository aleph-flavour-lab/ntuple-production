#ifndef ALEPHPHIKK_H
#define ALEPHPHIKK_H

/*
  phi(1020) -> K+ K- reconstruction from the FULL baseline-selected track list;
  no dE/dx quantity enters any selection. Input = flipD0_copy'ed trackstates
  (physical charge = +sign(omega)); lengths cm, momenta GeV.
*/

#include <algorithm>
#include <cmath>

#include <ROOT/RVec.hxx>
#include "TVector3.h"

#include "edm4hep/TrackState.h"
#include "FCCAnalyses/VertexingUtils.h"
#include "FCCAnalyses/VertexFitterSimple.h"

#include "aleph_units.h"
#include "analyzer_trkaux.h"
#include "analyzer_v0new.h"

namespace FCCAnalyses {
namespace AlephPhiKK {

using ROOT::VecOps::RVec;

// shared track auxiliaries (see analyzer_trkaux.h)
using AlephTrkAux::perigeeMomentum;
using AlephTrkAux::vertexDistSig;
using AlephTrkAux::apVars;
using AlephTrkAux::fitTracksCm;
using AlephTrkAux::TrkBlock;
using AlephTrkAux::TrkAux;
using AlephTrkAux::pushTrk;

constexpr double M_K   = AlephMasses::kK;    // charged kaon
constexpr double M_PHI = AlephMasses::kPhi;

// daughter momentum and energy in the phi rest frame
inline double pstarKK() {
  return std::sqrt(0.25 * M_PHI * M_PHI - M_K * M_K);  // 0.1269181 GeV
}
constexpr double estarKK() { return 0.5 * M_PHI; }

// Storage defaults, deliberately loose: the sample is re-cut offline.
constexpr double M_LO = 0.98, M_HI = 1.10;   // stored KK mass window [GeV]
constexpr double PRE_MARGIN = 0.02;          // pre-fit window margin [GeV]
constexpr double CHI2_CUT = 25.;             // loose sanity only (ndf = 1)
constexpr double DPV_FID = 50.;              // |vtx-PV| storage fiducial [cm]
// applied to the PERIGEE momentum, while the stored trk p is the at-vertex one
constexpr double P_MIN_DEF = 0.3;            // per-track |p| floor [GeV]

// Working-point flags, evaluated on the stored quantities; charge-blind.
constexpr double WP_DM = 0.012, TIGHT_DM = 0.005;  // |m - m_phi| [GeV]
constexpr double WP_PDAU = 1.0;                    // daughter |p| [GeV]
constexpr double WP_DPV = 1.0;                     // |vtx-PV| [cm]
constexpr double TIGHT_SIGD0 = 0.01;               // daughter sigma(d0) [cm]

// Armenteros-Podolanski band variable, 1 on the exact phi -> K+K- ellipse.
inline double phiBandEll(double alpha, double qt, double pmag) {
  const double ps = pstarKK();
  double beta = pmag / std::sqrt(pmag * pmag + M_PHI * M_PHI);
  if (beta <= 0.) return -1.;
  double amax = ps / (beta * estarKK());
  return std::sqrt(alpha * alpha / (amax * amax) + qt * qt / (ps * ps));
}

// One entry per accepted track pair; a track may appear in several candidates.
struct PhiKKCands {
  RVec<float> invM;       // KK invariant mass at the fitted vertex [GeV]
  RVec<float> p, px, py, pz;  // pair momentum at the fitted vertex [GeV]
  RVec<float> alpha, qt;  // Armenteros-Podolanski (qt in GeV)
  RVec<float> bandEll;    // equal-mass AP band variable (1 = exact locus)
  RVec<float> chi2;       // vertex fit chi2, normalised (ndf = 1)
  RVec<float> vx, vy, vz; // fitted vertex position [cm]
  RVec<float> dpv;        // |vertex - PV| [cm]
  RVec<float> dpvSig;     // 3D significance of the same, -1 if undefined
  RVec<int>   same_sign;  // 1 = same-charge pair (combinatorial control)
  // per-daughter blocks; trk1 is the higher-|p| track of the pair
  TrkBlock    trk1, trk2;
  RVec<int>   wp, tight;  // working-point flags, charge-blind
};

// The finder; every selection value is a constant defined above.
// veto_orig = original-Tracks indices excluded from the pairing pool.
inline PhiKKCands findPhiKK(
    const RVec<edm4hep::TrackState>& tracks,
    const RVec<int>& orig_idx,
    const RVec<int>& nvdet,
    const RVec<int>& nitc,
    const RVec<float>& chi2ndf,
    const RVec<int>& isprim,
    const VertexingUtils::FCCAnalysesVertex& PV,
    double solenoidBz,
    const RVec<int>& veto_orig,
    double bsx, double bsy, double bsz) {

  PhiKKCands out;
  const int nTr = tracks.size();
  if (nTr < 2) return out;

  const TrkAux aux{&orig_idx, &nvdet, &nitc, &chi2ndf, &isprim, nullptr, bsx, bsy, bsz};
  auto oidx = [&](int k) {
    return (k >= 0 && k < (int)orig_idx.size()) ? orig_idx[k] : -1;
  };

  // momentum prefilter; eK = kaon-hypothesis energy at the perigee
  const std::vector<char> vetoed = AlephTrkAux::memberMask(veto_orig);
  std::vector<int> good;
  good.reserve(nTr);
  std::vector<TVector3> pper(nTr);
  std::vector<double> eK(nTr, 0.);
  for (int i = 0; i < nTr; ++i) {
    const auto& t = tracks[i];
    if (AlephTrkAux::inMask(vetoed, oidx(i))) continue;
    pper[i] = perigeeMomentum(t, solenoidBz);
    if (pper[i].Mag() < P_MIN_DEF) continue;
    eK[i] = std::sqrt(pper[i].Mag2() + M_K * M_K);
    good.push_back(i);
  }
  if (good.size() < 2) return out;

  TVector3 pv(PV.vertex.position[0], PV.vertex.position[1], PV.vertex.position[2]);
  const double pre_lo = M_LO - PRE_MARGIN, pre_hi = M_HI + PRE_MARGIN;

  RVec<edm4hep::TrackState> tr_pair(2);
  for (size_t a = 0; a + 1 < good.size(); ++a) {
    const int i = good[a];
    for (size_t b = a + 1; b < good.size(); ++b) {
      const int j = good[b];
      const bool ss = (tracks[i].omega * tracks[j].omega > 0);

      // pre-fit KK mass from the perigee momenta
      TVector3 p_pre = pper[i] + pper[j];
      double e_pre = eK[i] + eK[j];
      double m_pre = std::sqrt(std::max(0., e_pre * e_pre - p_pre.Mag2()));
      if (m_pre < pre_lo || m_pre > pre_hi) continue;

      tr_pair[0] = tracks[i];
      tr_pair[1] = tracks[j];
      auto v = fitTracksCm(tr_pair, tracks, solenoidBz);
      if (v.updated_track_momentum_at_vertex.size() != 2) continue;

      double chi2 = v.vertex.chi2;
      if (!(chi2 == chi2)) continue;
      if (chi2 >= CHI2_CUT) continue;

      TVector3 pa = v.updated_track_momentum_at_vertex[0];
      TVector3 pb = v.updated_track_momentum_at_vertex[1];
      double m = AlephV0New::invMass(pa, M_K, pb, M_K);
      if (m < M_LO || m > M_HI) continue;

      TVector3 p = pa + pb;
      double pmag = p.Mag();
      if (!(pmag > 0.)) continue;
      // physical charge = +sign(omega); a same-sign pair follows track order
      int qi = (tracks[i].omega > 0) ? 1 : -1;
      double alpha, qt;
      apVars(pa, pb, ss ? 1 : qi, alpha, qt);
      double ell = phiBandEll(alpha, qt, pmag);

      TVector3 x(v.vertex.position[0], v.vertex.position[1], v.vertex.position[2]);
      TVector3 d = x - pv;
      double dis = d.Mag();
      if (dis > DPV_FID) continue;
      float dsig = vertexDistSig(d, v.vertex.covMatrix, PV.vertex.covMatrix);

      // daughter 1 = higher-|p| track
      int i1 = i, i2 = j;
      TVector3 q1 = pa, q2 = pb;
      if (pb.Mag() > pa.Mag()) { i1 = j; i2 = i; q1 = pb; q2 = pa; }

      out.invM.push_back(m);
      out.p.push_back(pmag);
      out.px.push_back(p.X()); out.py.push_back(p.Y()); out.pz.push_back(p.Z());
      out.alpha.push_back(alpha);
      out.qt.push_back(qt);
      out.bandEll.push_back(ell);
      out.chi2.push_back(chi2);
      out.vx.push_back(x.X()); out.vy.push_back(x.Y()); out.vz.push_back(x.Z());
      out.dpv.push_back(dis);
      out.dpvSig.push_back(dsig);
      out.same_sign.push_back(ss ? 1 : 0);
      pushTrk(out.trk1, tracks, i1, q1, aux);
      pushTrk(out.trk2, tracks, i2, q2, aux);

      const double dm = std::abs(m - M_PHI);
      const bool pdau = (q1.Mag() > WP_PDAU && q2.Mag() > WP_PDAU);
      const bool prompt = (dis < WP_DPV);
      const double sd1 = out.trk1.sigd0.back(), sd2 = out.trk2.sigd0.back();
      out.wp.push_back((dm < WP_DM && pdau && prompt) ? 1 : 0);
      out.tight.push_back((dm < TIGHT_DM && pdau && prompt &&
                           sd1 > 0. && sd1 < TIGHT_SIGD0 &&
                           sd2 > 0. && sd2 < TIGHT_SIGD0) ? 1 : 0);
    }
  }
  return out;
}

} // namespace AlephPhiKK
} // namespace FCCAnalyses

#endif
