#ifndef ALEPHV0NEW_H
#define ALEPHV0NEW_H

/*
  Standalone improved V0 reconstruction; output type
  VertexingUtils::FCCAnalysesV0. Input = flipD0_copy'ed trackstates (params AND
  covariance ALEPH->LCIO transformed); all fit-chain positions are in cm.
*/

#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <vector>

#include <ROOT/RVec.hxx>
#include "TVector3.h"

#include "aleph_units.h"
#include "analyzer_trkaux.h"
#include "dedx_valid.h"
#include "edm4hep/TrackData.h"
#include "edm4hep/TrackState.h"
#include "FCCAnalyses/VertexingUtils.h"
#include "FCCAnalyses/VertexFitterSimple.h"

namespace FCCAnalyses {
namespace AlephV0New {

using ROOT::VecOps::RVec;
using AlephTrkAux::apVars;
using AlephTrkAux::fitTracksCm;

constexpr double m_pi_ = AlephMasses::kPiCh;
constexpr double m_p_  = AlephMasses::kProton;
constexpr double MKS   = AlephMasses::kKs;
constexpr double MLAM  = AlephMasses::kLambda;

// ---------------------------------------------------------------------------
// Cut package: single named source. TIGHT = adopted package (findV0s defaults,
// stored per candidate as V0Collection::tight), LOOSE = ML-training superset
// tier. Mass windows, chi2 and displacement window are COMMON to both tiers.
// ---------------------------------------------------------------------------
constexpr double KS_M_LO = 0.40, KS_M_HI = 0.60;
constexpr double LAM_M_LO = 1.08, LAM_M_HI = 1.20;
constexpr double DIS_LO = 0.1, DIS_HI = 150.;
constexpr double CHI2_CUT = 10.;
constexpr double TIGHT_COS_KS_LOWP = 0.999, TIGHT_COS_KS_MIDP = 0.9995,
                 TIGHT_COS_KS_HIGHP = 0.9999;
// Lambda tight pointing, p-tiered in cos; deliberately NOT a mirror of the
// Ks ladder.
constexpr double TIGHT_COS_LAM_LOWP = 0.99995, TIGHT_COS_LAM_MIDP = 0.9999,
                 TIGHT_COS_LAM_HIGHP = 0.9999;
constexpr double TIGHT_QT_MIN_LAM = 0.04;
constexpr double AP_BAND_KS = 0.05, AP_LAM_LO = 0.10;
constexpr double TIGHT_NSIG_KS_LOWP = 3., TIGHT_NSIG_KS_HIGHP = 4.;
constexpr double LOOSE_COS_POINT = 0.999;
constexpr double LOOSE_QT_MIN_LAM = 0.02;
constexpr double LOOSE_NSIG_KS = 6.;
constexpr double LOOSE_LAM_BAND_LO = 0.20, LOOSE_LAM_BAND_HI = 0.40;
// LOOSE Lambda band-distance acceptance: the stored loose tier keeps
// |ell-1| / thr(p) below this fraction of the ramp half-width (see
// lamBandThrLoose).
constexpr double LOOSE_LAM_BAND_FRAC = 0.8;
constexpr double LAM_P_LO = 8., LAM_P_HI = 20.;
// Lambda AP-band ellipse resolution sigma_ell(p), quadrature model;
// "nsig" is in units of this 68% width.
constexpr double SIG_ELL_LAM_A = 0.01622, SIG_ELL_LAM_B = 0.0033748,
                 SIG_ELL_LAM_C = 0.00015544;
// Ks AP-band ellipse resolution, linear model — single source for the
// tight/loose band cuts AND the bandSig pull.
constexpr double SIG_ELL_KS_A = 0.007, SIG_ELL_KS_B = 0.0015;
constexpr double TIGHT_LAM_NSIG = 3.;
// Mass resolution sigma_m(p) [GeV], quadrature model — used by candMassSig.
constexpr double SIG_M_KS_A = 2.658e-3, SIG_M_KS_B = 0.5214e-3,
                 SIG_M_KS_C = 0.01418e-3;
constexpr double SIG_M_LAM_A = 1.045e-3, SIG_M_LAM_B = 0.2357e-3,
                 SIG_M_LAM_C = 0.005511e-3;

// Candidate-momentum [GeV] boundaries of the threshold tiers: the pointing cut
// steps at kPointTierLoP and kPointTierHiP, the AP band width at kBandTierP.
constexpr double kPointTierLoP = 2., kPointTierHiP = 4.;
constexpr double kBandTierP = 15.;

// Shared per-hypothesis acceptance helpers, used by findV0s for both tiers.
inline double ksPointThr(double pmag, double lowp, double midp, double highp) {
  return (pmag < kPointTierLoP) ? lowp : (pmag < kPointTierHiP) ? midp : highp;
}
// Lambda tight pointing tiers — same tier boundaries as Ks, different low-p
// value.
inline double lamPointThr(double pmag, double lowp = TIGHT_COS_LAM_LOWP,
                          double midp = TIGHT_COS_LAM_MIDP,
                          double highp = TIGHT_COS_LAM_HIGHP) {
  return (pmag < kPointTierLoP) ? lowp : (pmag < kPointTierHiP) ? midp : highp;
}
inline double ksBandEll(double alpha, double qt, double pmag) {
  // frozen tuned values of the Ks AP band; do not re-derive
  const double PSTAR_K = 0.20582, ESTAR_K = 0.248806;
  double beta = pmag / std::sqrt(pmag * pmag + MKS * MKS);
  double amax = PSTAR_K / (beta * ESTAR_K);
  return std::sqrt(std::pow(alpha / amax, 2) + std::pow(qt / PSTAR_K, 2));
}
inline double sigmaEllKs(double pmag) {
  return SIG_ELL_KS_A + SIG_ELL_KS_B * pmag;
}
inline double ksBandThr(double pmag, double floor_, double nsig_lo, double nsig_hi) {
  // resolution-scaled width; floor_ acts as the low-p floor
  double nsig = (pmag < kBandTierP) ? nsig_lo : nsig_hi;
  return std::max(floor_, nsig * sigmaEllKs(pmag));
}
inline double lamBandEll(double alpha, double qt, double pmag) {
  // frozen tuned values of the Lambda AP band; do not re-derive
  const double PSTAR_L = 0.1005, ALPHA0_L = 0.69157;
  double beta = pmag / std::sqrt(pmag * pmag + MLAM * MLAM);
  double amp = 2. * PSTAR_L / (beta * MLAM);
  return std::sqrt(std::pow((std::abs(alpha) - ALPHA0_L) / amp, 2) +
                   std::pow(qt / PSTAR_L, 2));
}
inline double lamBandThr(double pmag, double lo, double hi) {
  return (pmag < LAM_P_LO) ? lo : (pmag < LAM_P_HI) ? lo + (hi - lo) * (pmag - LAM_P_LO) / (LAM_P_HI - LAM_P_LO) : hi;
}
inline double sigmaEllLam(double pmag) {
  return std::sqrt(SIG_ELL_LAM_A * SIG_ELL_LAM_A +
                   std::pow(SIG_ELL_LAM_B * pmag, 2) +
                   std::pow(SIG_ELL_LAM_C * pmag * pmag, 2));
}
// TIGHT Lambda AP band: resolution-scaled, floored at floor_ and capped at the
// NOMINAL loose ramp edge (so the tight package is config-independent). The
// tight-inside-loose invariant is enforced on the loose side.
inline double lamBandThrTight(double pmag, double floor_ = AP_LAM_LO,
                              double nsig = TIGHT_LAM_NSIG) {
  return std::min(std::max(floor_, nsig * sigmaEllLam(pmag)),
                  lamBandThr(pmag, LOOSE_LAM_BAND_LO, LOOSE_LAM_BAND_HI));
}

// LOOSE Lambda AP band: BAND-DISTANCE convention — acceptance is a fixed
// fraction of the ramp half-width (|ell-1| < LOOSE_LAM_BAND_FRAC * thr(p)),
// floored at the tight threshold (tight_floor mirrors the caller's tight-clause
// floor) so the loose tier stays a superset of the tight one at every p.
inline double lamBandThrLoose(double pmag, double lo, double hi,
                              double tight_floor = AP_LAM_LO) {
  return std::max(LOOSE_LAM_BAND_FRAC * lamBandThr(pmag, lo, hi),
                  lamBandThrTight(pmag, tight_floor));
}

// TIGHT (adopted) package for ONE hypothesis: mass window, p-tiered pointing
// and AP band, plus the qT veto for Lambda. Single source for the finder tier.
inline bool ksTight(double m, double cp, double pmag, double alpha, double qt) {
  bool ok = (m > KS_M_LO && m < KS_M_HI) &&
            cp > ksPointThr(pmag, TIGHT_COS_KS_LOWP, TIGHT_COS_KS_MIDP,
                            TIGHT_COS_KS_HIGHP);
  if (ok)
    ok = std::abs(ksBandEll(alpha, qt, pmag) - 1.) <
         ksBandThr(pmag, AP_BAND_KS, TIGHT_NSIG_KS_LOWP, TIGHT_NSIG_KS_HIGHP);
  return ok;
}

inline bool lamTight(double m, double cp, double pmag, double alpha, double qt) {
  bool ok = (m > LAM_M_LO && m < LAM_M_HI) && cp > lamPointThr(pmag) &&
            qt > TIGHT_QT_MIN_LAM;
  if (ok)
    ok = std::abs(lamBandEll(alpha, qt, pmag) - 1.) < lamBandThrTight(pmag);
  return ok;
}

// momenta of the two tracks at the fitted vertex, already rescaled to the true
// GeV scale inside findV0s
inline void pairMomenta(const VertexingUtils::FCCAnalysesVertex& v,
                        TVector3& p1, TVector3& p2) {
  p1 = v.updated_track_momentum_at_vertex[0];
  p2 = v.updated_track_momentum_at_vertex[1];
}

inline double invMass(const TVector3& p1, double m1, const TVector3& p2, double m2) {
  double e1 = std::sqrt(p1.Mag2() + m1 * m1);
  double e2 = std::sqrt(p2.Mag2() + m2 * m2);
  TVector3 p = p1 + p2;
  double e = e1 + e2;
  return std::sqrt(std::max(0., e * e - p.Mag2()));
}

// ---------------------------------------------------------------------------
// Standard selection of ONE fitted pair, from quantities derived from the fit.
// tier: 0 rejected, 1 loose, 2 tight; pdg/m = booked hypothesis when tier > 0.
// ---------------------------------------------------------------------------
struct V0Sel { int tier; int pdg; double m; };

inline V0Sel evalV0Selection(double chi2, double dis, double cp, double pmag,
                             double alpha, double qt, double mks, double mlam) {
  V0Sel s{0, 310, mks};
  if (chi2 >= CHI2_CUT || !(chi2 == chi2)) return s;
  if (dis < DIS_LO || dis > DIS_HI) return s;
  if (pmag <= 0) return s;

  // TIGHT (adopted) package first; arbitration among the tight-passing
  // hypotheses only, so the tight subset is EXACTLY what the module would
  // output with the loose tier switched off.
  bool okKs = ksTight(mks, cp, pmag, alpha, qt);
  bool okLam = lamTight(mlam, cp, pmag, alpha, qt);
  bool tight = okKs || okLam;
  if (!tight) {
    // LOOSE training tier: flat pointing, widened AP bands, relaxed
    // Lambda qT veto; windows/chi2/displacement common.
    bool inWinKs = (mks > KS_M_LO && mks < KS_M_HI);
    bool inWinLam = (mlam > LAM_M_LO && mlam < LAM_M_HI);
    okKs = inWinKs && cp > LOOSE_COS_POINT;
    if (okKs)
      okKs = std::abs(ksBandEll(alpha, qt, pmag) - 1.) <
             ksBandThr(pmag, AP_BAND_KS, LOOSE_NSIG_KS, LOOSE_NSIG_KS);
    okLam = inWinLam && cp > LOOSE_COS_POINT && qt > LOOSE_QT_MIN_LAM;
    if (okLam)
      okLam = std::abs(lamBandEll(alpha, qt, pmag) - 1.) <
              lamBandThrLoose(pmag, LOOSE_LAM_BAND_LO, LOOSE_LAM_BAND_HI, AP_LAM_LO);
    if (!okKs && !okLam) return s;
  }
  double dks = std::abs(mks - MKS) / (0.5 * (KS_M_HI - KS_M_LO));
  double dlam = std::abs(mlam - MLAM) / (0.5 * (LAM_M_HI - LAM_M_LO));
  if (okKs && (!okLam || dks <= dlam)) { s.pdg = 310;  s.m = mks; }
  else                                 { s.pdg = 3122; s.m = mlam; }
  s.tier = tight ? 2 : 1;
  return s;
}

// ---------------------------------------------------------------------------
// The finder. np_tracks = flipD0_copy'ed non-primary trackstates, PV = fitted
// primary vertex (positions in cm). TWO-TIER: only tight-failing pairs enter the
// LOOSE tier, and tight candidates claim tracks first, so the tight flag selects
// exactly the tight-only output. Returns candidates in claim order (tight block
// first, chi2 ascending within a tier); pdgAbs = best hypothesis (310 or 3122),
// invM its mass, tight = 1 for the tight tier.
// ---------------------------------------------------------------------------

// The module's candidate collection: the FCCAnalyses V0 collection (vertices,
// hypothesis, mass), so every helper taking one accepts it, plus the tier of
// each candidate. tightV0s / getKs / getLambda copy sub-collections of it.
struct V0Collection : VertexingUtils::FCCAnalysesV0 {
  RVec<int> tight;   // 1 = passed the tight package, 0 = loose tier
};

inline V0Collection findV0s(
    const RVec<edm4hep::TrackState>& np_tracks,
    const VertexingUtils::FCCAnalysesVertex& PV,
    double solenoidBz) {

  V0Collection result;
  const int nTr = np_tracks.size();
  if (nTr < 2) return result;

  TVector3 pv(PV.vertex.position[0], PV.vertex.position[1], PV.vertex.position[2]);

  struct Cand {
    VertexingUtils::FCCAnalysesVertex vtx;
    int i, j;
    int pdg;       // best hypothesis
    double m;      // mass under best hypothesis
    double chi2;
    bool tight;    // passed the tight (adopted) package
  };
  std::vector<Cand> cands;

  RVec<edm4hep::TrackState> tr_pair(2);
  for (int i = 0; i < nTr - 1; ++i) {
    tr_pair[0] = np_tracks[i];
    for (int j = i + 1; j < nTr; ++j) {
      if (np_tracks[i].omega * np_tracks[j].omega > 0) continue; // same charge
      tr_pair[1] = np_tracks[j];

      auto v = fitTracksCm(tr_pair, np_tracks, solenoidBz);
      if (v.updated_track_momentum_at_vertex.size() != 2) continue;
      double chi2 = v.vertex.chi2; // normalised, ndf=1

      // displacement (cm) + pointing
      TVector3 x(v.vertex.position[0], v.vertex.position[1], v.vertex.position[2]);
      TVector3 d = x - pv;
      double dis = d.Mag();

      TVector3 p1, p2;
      pairMomenta(v, p1, p2);
      TVector3 p = p1 + p2;
      double pmag = p.Mag();
      double cp = (dis > 0. && pmag > 0.) ? d.Dot(p) / (dis * pmag) : -2.;
      // physical charge = +sign(omega) for the flipD0_copy'ed collection (raw
      // ALEPH omega carries -charge, the flip restores +charge)
      double q1 = (np_tracks[i].omega > 0) ? 1. : -1.;
      double alpha = AlephTrkAux::kApUndef, qt = AlephTrkAux::kApUndef;
      if (pmag > 0.) apVars(p1, p2, q1, alpha, qt);

      // hypothesis masses: Ks(pipi), Lambda(p pi) with proton = higher-|p| track
      // (in a Lambda decay the baryon carries most of the momentum)
      double mks = invMass(p1, m_pi_, p2, m_pi_);
      double mlam = (p1.Mag() > p2.Mag()) ? invMass(p1, m_p_, p2, m_pi_)
                                          : invMass(p1, m_pi_, p2, m_p_);

      V0Sel sel = evalV0Selection(chi2, dis, cp, pmag, alpha, qt, mks, mlam);
      if (sel.tier == 0) continue;

      cands.push_back({v, i, j, sel.pdg, sel.m, chi2, sel.tier == 2});
    }
  }

  // quality-ranked global claiming: tight candidates claim first (preserving
  // the tight-only output), then loose; best chi2 first within a tier
  std::vector<size_t> order(cands.size());
  std::iota(order.begin(), order.end(), 0);
  std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) {
    if (cands[a].tight != cands[b].tight) return cands[a].tight;
    return cands[a].chi2 < cands[b].chi2;
  });

  std::vector<bool> used(nTr, false);
  for (size_t k : order) {
    const Cand& c = cands[k];
    if (used[c.i] || used[c.j]) continue;
    used[c.i] = true;
    used[c.j] = true;
    result.vtx.push_back(c.vtx);
    result.pdgAbs.push_back(c.pdg);
    result.invM.push_back(c.m);
    result.tight.push_back(c.tight ? 1 : 0);
  }
  return result;
}

// Sub-collections, in the stored order: the tight tier, or the candidates
// booked on one hypothesis (310 = Ks, 3122 = Lambda or anti-Lambda). Copies of
// the same fitted objects, no refit: getKs(V0sNewTight_event) is the tight Ks,
// getKs(V0sNew_event) every stored Ks candidate (loose tier, tight included).
inline V0Collection selectV0(const V0Collection& v0s, const RVec<int>& keep) {
  if (keep.size() != v0s.vtx.size())
    throw std::runtime_error("AlephV0New::selectV0: mask size != candidate count");
  V0Collection out;
  for (size_t c = 0; c < v0s.vtx.size(); ++c) {
    if (!keep[c]) continue;
    out.vtx.push_back(v0s.vtx[c]);
    out.pdgAbs.push_back(v0s.pdgAbs[c]);
    out.invM.push_back(v0s.invM[c]);
    out.tight.push_back(v0s.tight[c]);
  }
  return out;
}
inline V0Collection tightV0s(const V0Collection& v0s) {
  return selectV0(v0s, v0s.tight);
}
inline V0Collection getKs(const V0Collection& v0s) {
  return selectV0(v0s, v0s.pdgAbs == 310);
}
inline V0Collection getLambda(const V0Collection& v0s) {
  return selectV0(v0s, v0s.pdgAbs == 3122);
}

// ---------------------------------------------------------------------------
// Truth-free per-candidate diagnostics (work on data; reco_ind is filled by
// this module, momenta are already at the true GeV scale).
// ---------------------------------------------------------------------------

inline RVec<float> candAlpha(const VertexingUtils::FCCAnalysesV0& v0s,
                             const RVec<edm4hep::TrackState>& secondaries) {
  RVec<float> out;
  for (const auto& v : v0s.vtx) {
    if (v.reco_ind.size() < 2 || v.updated_track_momentum_at_vertex.size() < 2 ||
        v.reco_ind[0] < 0 || v.reco_ind[0] >= (int)secondaries.size()) {
      out.push_back(AlephTrkAux::kApUndef);
      continue;
    }
    double q1 = (secondaries[v.reco_ind[0]].omega > 0) ? 1. : -1.;
    double alpha, qt;
    apVars(v.updated_track_momentum_at_vertex[0],
           v.updated_track_momentum_at_vertex[1], q1, alpha, qt,
           AlephTrkAux::kApUndef);
    out.push_back(alpha);
  }
  return out;
}

// ML-input pulls of the BOOKED hypothesis, in resolution units. Both SIGNED;
// -999 = undefined candidate. bandSig = (bandEll - 1) / sigma_ell(p),
// massSig = (invM - m_hyp) / sigma_m(p). Other cut variables are stored raw.
inline RVec<float> candBandSig(const VertexingUtils::FCCAnalysesV0& v0s,
                               const RVec<edm4hep::TrackState>& secondaries) {
  RVec<float> out;
  for (size_t c = 0; c < v0s.vtx.size(); ++c) {
    const auto& v = v0s.vtx[c];
    if (v.reco_ind.size() < 2 || v.updated_track_momentum_at_vertex.size() < 2 ||
        v.reco_ind[0] < 0 || v.reco_ind[0] >= (int)secondaries.size()) {
      out.push_back(-999.);
      continue;
    }
    TVector3 p1 = v.updated_track_momentum_at_vertex[0];
    TVector3 p2 = v.updated_track_momentum_at_vertex[1];
    TVector3 p = p1 + p2;
    double pmag = p.Mag();
    if (pmag <= 0) { out.push_back(-999.); continue; }
    double q1 = (secondaries[v.reco_ind[0]].omega > 0) ? 1. : -1.;
    double alpha, qt;
    apVars(p1, p2, q1, alpha, qt);
    if (v0s.pdgAbs[c] == 310)
      out.push_back((ksBandEll(alpha, qt, pmag) - 1.) / sigmaEllKs(pmag));
    else
      out.push_back((lamBandEll(alpha, qt, pmag) - 1.) / sigmaEllLam(pmag));
  }
  return out;
}

inline RVec<float> candMassSig(const VertexingUtils::FCCAnalysesV0& v0s) {
  RVec<float> out;
  for (size_t c = 0; c < v0s.vtx.size(); ++c) {
    const auto& v = v0s.vtx[c];
    if (v.updated_track_momentum_at_vertex.size() < 2) {
      out.push_back(-999.);
      continue;
    }
    TVector3 p = v.updated_track_momentum_at_vertex[0] +
                 v.updated_track_momentum_at_vertex[1];
    double pmag = p.Mag(), p2 = pmag * pmag;
    bool isKs = (v0s.pdgAbs[c] == 310);
    double a = isKs ? SIG_M_KS_A : SIG_M_LAM_A;
    double b = isKs ? SIG_M_KS_B : SIG_M_LAM_B;
    double cc = isKs ? SIG_M_KS_C : SIG_M_LAM_C;
    double sig = std::sqrt(a * a + b * b * p2 + cc * cc * p2 * p2);
    out.push_back((v0s.invM[c] - (isKs ? MKS : MLAM)) / sig);
  }
  return out;
}

// Pointing significance: chi2-like significance of the displacement component
// PERPENDICULAR to the candidate momentum (all in cm). d = candidate vertex -
// reference vertex, p = candidate momentum, cV/cR = packed lower-triangular
// position covariances (xx,yx,yy,zx,zy,zz); the reference vertex is the PV.
// Returns -1 for degenerate or singular geometry.
template <typename CovV, typename CovR>
inline float pointSigTransverse(const TVector3& d, const TVector3& p,
                                const CovV& cV, const CovR& cR) {
  if (p.Mag() <= 0 || d.Mag() <= 0) return -1.;
  TVector3 ph = p.Unit();
  TVector3 u1 = ph.Orthogonal().Unit();
  TVector3 u2 = ph.Cross(u1);
  double C[3][3];
  AlephTrkAux::sumCovPacked(cV, cR, C);
  auto quad = [&](const TVector3& a, const TVector3& b) {
    double s = 0.;
    double av[3] = {a.X(), a.Y(), a.Z()}, bv[3] = {b.X(), b.Y(), b.Z()};
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) s += av[i] * C[i][j] * bv[j];
    return s;
  };
  double c11 = quad(u1, u1), c22 = quad(u2, u2), c12 = quad(u1, u2);
  double det = c11 * c22 - c12 * c12;
  if (det <= 0. || c11 <= 0. || c22 <= 0.) return -1.;
  double d1 = d.Dot(u1), d2 = d.Dot(u2);
  double sig2 = (d1 * (c22 * d1 - c12 * d2) + d2 * (c11 * d2 - c12 * d1)) / det;
  return sig2 > 0. ? std::sqrt(sig2) : 0.;
}

inline RVec<float> candPointSig(const VertexingUtils::FCCAnalysesV0& v0s,
                                const VertexingUtils::FCCAnalysesVertex& PV) {
  RVec<float> out;
  TVector3 pv(PV.vertex.position[0], PV.vertex.position[1], PV.vertex.position[2]);
  for (const auto& v : v0s.vtx) {
    TVector3 x(v.vertex.position[0], v.vertex.position[1], v.vertex.position[2]);
    TVector3 p(0., 0., 0.);
    for (const auto& tp : v.updated_track_momentum_at_vertex) p += tp;
    out.push_back(pointSigTransverse(x - pv, p, v.vertex.covMatrix,
                                     PV.vertex.covMatrix));
  }
  return out;
}

inline RVec<float> candQt(const VertexingUtils::FCCAnalysesV0& v0s) {
  RVec<float> out;
  for (const auto& v : v0s.vtx) {
    TVector3 pa = v.updated_track_momentum_at_vertex[0];
    TVector3 p = pa + v.updated_track_momentum_at_vertex[1];
    out.push_back(pa.Cross(p.Unit()).Mag());
  }
  return out;
}

// ---------------------------------------------------------------------------
// vertex-fit covariance exposure + per-daughter joins
// ---------------------------------------------------------------------------

// Vertex-fit covariance component ic of every candidate, in cm^2 (packed lower
// triangle: 0=xx 1=yx 2=yy 3=zx 4=zy 5=zz).
inline RVec<float> candCovComp(const VertexingUtils::FCCAnalysesV0& v0s,
                               int ic) {
  RVec<float> out;
  for (const auto& v : v0s.vtx) out.push_back(v.vertex.covMatrix[ic]);
  return out;
}

// Slot of logical daughter k (0/1) in the parallel per-track arrays reco_ind
// and updated_track_momentum_at_vertex (same slot = same fitted track): leg 0 =
// higher-momentum daughter at the fitted vertex. Raw order on a tie or when the
// momenta are missing.
inline int legSlot(const VertexingUtils::FCCAnalysesVertex& v, int k) {
  if (v.updated_track_momentum_at_vertex.size() < 2) return k;
  const bool swap = v.updated_track_momentum_at_vertex[1].Mag() >
                    v.updated_track_momentum_at_vertex[0].Mag();
  return swap ? 1 - k : k;
}

// Daughter k (0/1) of every candidate as an ORIGINAL Tracks index: reco_ind
// (secondary space) walked through sec2orig, legs in momentum order (k = 0 is
// the higher-momentum one). -1 when unavailable; truth-free.
inline RVec<int> candDaughterOrigIdx(const VertexingUtils::FCCAnalysesV0& v0s,
                                     const RVec<int>& sec2orig, int k) {
  RVec<int> out;
  for (const auto& v : v0s.vtx) {
    int idx = -1;
    const int slot = legSlot(v, k);
    if (slot < (int)v.reco_ind.size()) {
      int s = v.reco_ind[slot];
      if (s >= 0 && s < (int)sec2orig.size()) idx = sec2orig[s];
    }
    out.push_back(idx);
  }
  return out;
}

// Measurement index of every original track index, or -1: the dE/dx join built
// ONCE per collection per event, in place of a scan per requested track. A
// track measured twice keeps the LAST measurement, as the particle-flow join
// does. That entry is then -1 when it fails the shared validity gate (value ==
// omega of the track = the failed-leg sentinel, or non-finite/non-positive
// value or error), so value and error branches share one lookup; gate = false
// accepts every linked measurement. Sized by the track count: callers index the
// result by original Tracks index.
inline RVec<int> dedxIndexByTrack(const RVec<float>& value,
                                  const RVec<float>& error,
                                  const RVec<int>& meas_track_idx,
                                  const RVec<edm4hep::TrackData>& tracks,
                                  const RVec<edm4hep::TrackState>& trackStates,
                                  bool gate = true) {
  RVec<int> out(tracks.size(), -1);
  const size_t nm = std::min({value.size(), error.size(),
                              meas_track_idx.size()});
  for (size_t j = 0; j < nm; ++j) {
    const int t = meas_track_idx[j];
    // out-of-range relation entry = corrupt input: fail loudly
    if (t < 0 || t >= (int)out.size())
      throw std::runtime_error(
          "AlephV0New::dedxIndexByTrack: dEdx->Track relation entry out of range");
    out[t] = (int)j;
  }
  if (!gate) return out;
  for (size_t t = 0; t < out.size(); ++t) {
    const int j = out[t];
    if (j < 0) continue;
    const size_t stateIdx = tracks[t].trackStates_begin;
    if (stateIdx >= trackStates.size())
      throw std::runtime_error(
          "AlephV0New::dedxIndexByTrack: Track->TrackState index out of range");
    if (!AlephDedx::dEdxValid(value[j], error[j], trackStates[stateIdx].omega))
      out[t] = -1;
  }
  return out;
}

// Per-track quantity by ORIGINAL track index through that join; -1 when the
// track has no valid measurement.
inline RVec<float> trackQuantityByIndex(const RVec<int>& want,
                                        const RVec<float>& values,
                                        const RVec<int>& meas_of_track) {
  RVec<float> out;
  for (int w : want) {
    float val = -1.f;
    if (w >= 0 && w < (int)meas_of_track.size()) {
      const int j = meas_of_track[w];
      if (j >= 0 && j < (int)values.size()) val = values[j];
    }
    out.push_back(val);
  }
  return out;
}

// ---------------------------------------------------------------------------
// Jet-relative and per-jet tagger inputs. Everything below is derived from the
// stored candidates, the primary vertex and the jets: no candidate is re-fitted
// and no new tuned value enters.
// ---------------------------------------------------------------------------

// Undefined value of the jet-relative float branches, and of the per-jet int
// branches that carry a flag rather than a count.
constexpr float TAG_UNDEF = AlephTrkAux::kUndef;
constexpr int TAG_UNDEF_INT = -1;

// Summed daughter momentum at the fitted vertex [GeV] - the momentum the jet
// assignment and every jet-relative quantity below use.
inline TVector3 candMomentum(const VertexingUtils::FCCAnalysesVertex& v) {
  TVector3 p(0., 0., 0.);
  for (const auto& tp : v.updated_track_momentum_at_vertex) p += tp;
  return p;
}

// Jet of every candidate, or -1 when it has none: closest dR between the
// candidate momentum and the jet axis, the first jet winning a tie. Reproduces
// the assignment that fills the per-jet mirror block (assign_V0s_to_jets in
// analyzer.h: same zero-momentum guard and dR seed), so the two are joinable.
inline RVec<int> candJetIdx(const VertexingUtils::FCCAnalysesV0& v0s,
                            const RVec<fastjet::PseudoJet>& jets) {
  RVec<int> out;
  for (const auto& v : v0s.vtx) {
    int best = -1;
    TVector3 p = candMomentum(v);
    if (jets.size() > 0 && !(p.Mag() < 1e-10)) {
      double minDR = 99.;
      best = 0;
      for (size_t j = 0; j < jets.size(); ++j) {
        double dR = p.DeltaR(TVector3(jets[j].px(), jets[j].py(), jets[j].pz()));
        if (dR < minDR) { minDR = dR; best = (int)j; }
      }
    }
    out.push_back(best);
  }
  return out;
}

// Jet-relative candidate kinematics, selected by `which`:
// 0 = |p| / E_jet, 1 = (p . jet direction) / E_jet, 2 = |p x jet direction|
// [GeV], 3 = dR between the candidate momentum and the jet axis. TAG_UNDEF when
// the candidate has no jet or the jet has no direction/energy.
inline RVec<float> candJetVar(const VertexingUtils::FCCAnalysesV0& v0s,
                              const RVec<fastjet::PseudoJet>& jets,
                              const RVec<int>& jetIdx, int which) {
  RVec<float> out;
  for (size_t c = 0; c < v0s.vtx.size(); ++c) {
    float val = TAG_UNDEF;
    const int j = (c < jetIdx.size()) ? jetIdx[c] : -1;
    if (j >= 0 && j < (int)jets.size()) {
      TVector3 ax(jets[j].px(), jets[j].py(), jets[j].pz());
      const double ej = jets[j].e();
      if (ax.Mag() > 0.) {
        TVector3 p = candMomentum(v0s.vtx[c]);
        TVector3 jh = ax.Unit();
        if (which == 0)      { if (ej > 0.) val = p.Mag() / ej; }
        else if (which == 1) { if (ej > 0.) val = p.Dot(jh) / ej; }
        else if (which == 2) { val = p.Cross(jh).Mag(); }
        else                 { val = p.DeltaR(ax); }
      }
    }
    out.push_back(val);
  }
  return out;
}

// Momentum rank of every candidate among the candidates of its own jet,
// 1 = leading; -1 for a candidate without a jet. An exact momentum tie is
// broken by candidate order.
inline RVec<int> candRankInJet(const VertexingUtils::FCCAnalysesV0& v0s,
                               const RVec<int>& jetIdx) {
  const size_t n = v0s.vtx.size();
  RVec<int> out(n, -1);
  std::vector<double> pm(n, 0.);
  for (size_t c = 0; c < n; ++c) pm[c] = candMomentum(v0s.vtx[c]).Mag();
  for (size_t c = 0; c < n; ++c) {
    if (c >= jetIdx.size() || jetIdx[c] < 0) continue;
    int rank = 1;
    for (size_t k = 0; k < n && k < jetIdx.size(); ++k) {
      if (k == c || jetIdx[k] != jetIdx[c]) continue;
      if (pm[k] > pm[c] || (pm[k] == pm[c] && k < c)) ++rank;
    }
    out[c] = rank;
  }
  return out;
}

// Transverse flight length of every candidate vertex from the reference vertex
// [cm]; the 3D one is candDxyz. TAG_UNDEF when the reference vertex has fewer
// than kPVMinTracks tracks, i.e. is the default vertex at the origin.
inline RVec<float> candLxy(const VertexingUtils::FCCAnalysesV0& v0s,
                           const VertexingUtils::FCCAnalysesVertex& PV) {
  RVec<float> out;
  if (PV.ntracks < AlephTrkAux::kPVMinTracks)
    return RVec<float>(v0s.vtx.size(), TAG_UNDEF);
  TVector3 pv(PV.vertex.position[0], PV.vertex.position[1], PV.vertex.position[2]);
  for (const auto& v : v0s.vtx) {
    TVector3 x(v.vertex.position[0], v.vertex.position[1], v.vertex.position[2]);
    out.push_back((x - pv).Perp());
  }
  return out;
}

// Flight-length significance L / sigma_L, with sigma_L the candidate and
// reference position covariances summed and projected on the flight direction:
// which 0 = transverse (xy), 1 = 3D. TAG_UNDEF for a vanishing flight length, a
// non-positive projected variance, or a reference vertex with fewer than
// kPVMinTracks tracks (the default vertex at the origin).
inline RVec<float> candFlightSig(const VertexingUtils::FCCAnalysesV0& v0s,
                                 const VertexingUtils::FCCAnalysesVertex& PV,
                                 int which) {
  RVec<float> out;
  if (PV.ntracks < AlephTrkAux::kPVMinTracks)
    return RVec<float>(v0s.vtx.size(), TAG_UNDEF);
  TVector3 pv(PV.vertex.position[0], PV.vertex.position[1], PV.vertex.position[2]);
  double C[3][3];
  for (const auto& v : v0s.vtx) {
    TVector3 x(v.vertex.position[0], v.vertex.position[1], v.vertex.position[2]);
    TVector3 d = x - pv;
    if (which == 0) d.SetZ(0.);
    const double L = d.Mag();
    if (!(L > 0.)) { out.push_back(TAG_UNDEF); continue; }
    AlephTrkAux::sumCovPacked(v.vertex.covMatrix, PV.vertex.covMatrix, C);
    TVector3 u = d.Unit();
    const double uv[3] = {u.X(), u.Y(), u.Z()};
    double s2 = 0.;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) s2 += uv[i] * C[i][j] * uv[j];
    out.push_back((s2 > 0. && std::isfinite(s2)) ? float(L / std::sqrt(s2))
                                                 : TAG_UNDEF);
  }
  return out;
}

// Baryon sign of every candidate: +1 Lambda, -1 anti-Lambda, 0 for a Ks or an
// undefined candidate. The proton is the higher-momentum leg, so the sign of
// alpha is the proton charge - the same convention candAlpha writes (kApUndef
// when undefined; |alpha| > 1 is physical when one leg points backwards). Any
// value within 1 of the sentinel counts as undefined.
inline RVec<int> candBaryon(const VertexingUtils::FCCAnalysesV0& v0s,
                            const RVec<edm4hep::TrackState>& secondaries) {
  RVec<int> out;
  const RVec<float> alpha = candAlpha(v0s, secondaries);
  for (size_t c = 0; c < v0s.pdgAbs.size(); ++c) {
    int b = 0;
    if (v0s.pdgAbs[c] == 3122 && c < alpha.size() &&
        alpha[c] > float(AlephTrkAux::kApUndef) + 1.f)
      b = (alpha[c] > 0.f) ? 1 : -1;
    out.push_back(b);
  }
  return out;
}

// Physical charge of daughter k (0/1) of every candidate, legs in momentum
// order: the secondary track states are flipD0_copy'ed, so the charge is
// +sign(omega). 0 when the daughter is unavailable.
inline RVec<int> candDaughterCharge(const VertexingUtils::FCCAnalysesV0& v0s,
                                    const RVec<edm4hep::TrackState>& secondaries,
                                    int k) {
  RVec<int> out;
  for (const auto& v : v0s.vtx) {
    int q = 0;
    const int slot = legSlot(v, k);
    if (slot < (int)v.reco_ind.size()) {
      const int s = v.reco_ind[slot];
      if (s >= 0 && s < (int)secondaries.size())
        q = (secondaries[s].omega > 0) ? 1 : -1;
    }
    out.push_back(q);
  }
  return out;
}

// Momentum magnitude of daughter k (0/1) at the fitted vertex [GeV], legs in
// momentum order (k = 0 is the higher-momentum one); TAG_UNDEF when the
// daughter is unavailable.
inline RVec<float> candDaughterP(const VertexingUtils::FCCAnalysesV0& v0s,
                                 int k) {
  RVec<float> out;
  for (const auto& v : v0s.vtx) {
    float p = TAG_UNDEF;
    const int slot = legSlot(v, k);
    if (slot < (int)v.updated_track_momentum_at_vertex.size())
      p = v.updated_track_momentum_at_vertex[slot].Mag();
    out.push_back(p);
  }
  return out;
}

// Daughters of every candidate (0, 1 or 2) whose original track is also a
// daughter of another stored candidate. Exclusive claiming makes this 0
// everywhere, so a non-zero entry flags a broken claim.
inline RVec<int> candNShared(const RVec<int>& d1, const RVec<int>& d2) {
  const size_t n = std::min(d1.size(), d2.size());
  std::map<int, int> nUse;
  for (size_t c = 0; c < n; ++c) {
    if (d1[c] >= 0) nUse[d1[c]]++;
    if (d2[c] >= 0) nUse[d2[c]]++;
  }
  RVec<int> out;
  for (size_t c = 0; c < n; ++c) {
    int s = 0;
    if (d1[c] >= 0 && nUse[d1[c]] > 1) ++s;
    if (d2[c] >= 0 && nUse[d2[c]] > 1) ++s;
    out.push_back(s);
  }
  return out;
}

// Candidates of each jet passing one species/tier combination. want_pdg = 310
// or 3122; want_baryon = +1/-1 to require that baryon sign, 0 = either;
// want_tight = 1 restricts to the tight tier, 0 counts every stored candidate.
inline RVec<int> jetCountV0(const RVec<int>& jetIdx,
                            const RVec<fastjet::PseudoJet>& jets,
                            const RVec<int>& pdg, const RVec<int>& tight,
                            const RVec<int>& baryon, int want_pdg,
                            int want_baryon, int want_tight) {
  RVec<int> out(jets.size(), 0);
  for (size_t c = 0; c < jetIdx.size(); ++c) {
    const int j = jetIdx[c];
    if (j < 0 || j >= (int)jets.size()) continue;
    if (c >= pdg.size() || pdg[c] != want_pdg) continue;
    if (want_tight && (c >= tight.size() || !tight[c])) continue;
    if (want_baryon != 0 && (c >= baryon.size() || baryon[c] != want_baryon))
      continue;
    ++out[j];
  }
  return out;
}

// Candidate index of the leading (highest-momentum) candidate of each jet, -1
// for a jet without candidates. Leading = rank 1 of candRankInJet.
inline RVec<int> jetLeadIdx(const RVec<int>& jetIdx,
                            const RVec<fastjet::PseudoJet>& jets,
                            const RVec<int>& rank) {
  RVec<int> out(jets.size(), -1);
  for (size_t c = 0; c < jetIdx.size() && c < rank.size(); ++c) {
    const int j = jetIdx[c];
    if (j >= 0 && j < (int)jets.size() && rank[c] == 1) out[j] = (int)c;
  }
  return out;
}

// Per-candidate quantity of each jet's leading candidate, `fill` for a jet
// without one; the default suits the code-valued branches, flags pass
// TAG_UNDEF_INT.
inline RVec<int> jetLeadInt(const RVec<int>& lead, const RVec<int>& values,
                            int fill = 0) {
  RVec<int> out;
  for (int c : lead)
    out.push_back((c >= 0 && c < (int)values.size()) ? values[c] : fill);
  return out;
}

inline RVec<float> jetLeadFloat(const RVec<int>& lead,
                                const RVec<float>& values, float fill) {
  RVec<float> out;
  for (int c : lead)
    out.push_back((c >= 0 && c < (int)values.size()) ? values[c] : fill);
  return out;
}

// Summed momentum fraction of the TIGHT candidates of each jet; 0 for a jet
// with none. Undefined fractions are left out of the sum.
inline RVec<float> jetSumZ(const RVec<int>& jetIdx,
                           const RVec<fastjet::PseudoJet>& jets,
                           const RVec<float>& z, const RVec<int>& tight) {
  RVec<float> out(jets.size(), 0.f);
  for (size_t c = 0; c < jetIdx.size() && c < z.size() && c < tight.size(); ++c) {
    const int j = jetIdx[c];
    if (j < 0 || j >= (int)jets.size() || !tight[c] || z[c] < 0.f) continue;
    out[j] += z[c];
  }
  return out;
}

// Distinct original tracks claimed by the TIGHT candidates of each jet.
inline RVec<int> jetNTrkClaimed(const RVec<int>& jetIdx,
                                const RVec<fastjet::PseudoJet>& jets,
                                const RVec<int>& tight, const RVec<int>& d1,
                                const RVec<int>& d2) {
  RVec<int> out(jets.size(), 0);
  std::vector<std::set<int>> claimed(jets.size());
  for (size_t c = 0; c < jetIdx.size(); ++c) {
    const int j = jetIdx[c];
    if (j < 0 || j >= (int)jets.size()) continue;
    if (c >= tight.size() || !tight[c]) continue;
    if (c < d1.size() && d1[c] >= 0) claimed[j].insert(d1[c]);
    if (c < d2.size() && d2[c] >= 0) claimed[j].insert(d2[c]);
  }
  for (size_t j = 0; j < jets.size(); ++j) out[j] = (int)claimed[j].size();
  return out;
}

// Per-candidate quantity regrouped per jet, in the order the per-jet mirror
// block is filled (candidate order within each jet), so the two are joinable.
inline RVec<RVec<int>> jetGatherInt(const RVec<int>& jetIdx,
                                    const RVec<fastjet::PseudoJet>& jets,
                                    const RVec<int>& values) {
  RVec<RVec<int>> out(jets.size());
  for (size_t c = 0; c < jetIdx.size(); ++c) {
    const int j = jetIdx[c];
    if (j < 0 || j >= (int)jets.size()) continue;
    out[j].push_back((c < values.size()) ? values[c] : -1);
  }
  return out;
}

// The candidate indices themselves, in the same per-jet order: the join key
// from the per-jet mirror block back to the event-level candidate list.
inline RVec<RVec<int>> jetGatherIdx(const RVec<int>& jetIdx,
                                    const RVec<fastjet::PseudoJet>& jets) {
  RVec<RVec<int>> out(jets.size());
  for (size_t c = 0; c < jetIdx.size(); ++c) {
    const int j = jetIdx[c];
    if (j >= 0 && j < (int)jets.size()) out[j].push_back((int)c);
  }
  return out;
}

} // namespace AlephV0New
} // namespace FCCAnalyses

#endif
