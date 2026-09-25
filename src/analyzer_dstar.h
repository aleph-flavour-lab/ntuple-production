#ifndef ALEPHDSTAR_H
#define ALEPHDSTAR_H

/*
  D*+ -> D0(K- pi+) pi+_slow reconstruction as a kinematically tagged kaon
  source; no dE/dx quantity enters any selection. Input = flipD0_copy'ed
  trackstates (physical charge = +sign(omega)); lengths cm, momenta GeV.
*/

#include <algorithm>
#include <cmath>
#include <vector>

#include <ROOT/RVec.hxx>
#include "TVector3.h"
#include "TLorentzVector.h"

#include "edm4hep/TrackState.h"
#include "FCCAnalyses/VertexingUtils.h"
#include "FCCAnalyses/VertexFitterSimple.h"

#include "aleph_units.h"
#include "analyzer_trkaux.h"
#include "analyzer_v0new.h"

namespace FCCAnalyses {
namespace AlephDstar {

using ROOT::VecOps::RVec;

// shared track auxiliaries (see analyzer_trkaux.h)
using AlephTrkAux::perigeeMomentum;
using AlephTrkAux::vertexDistSig;
using AlephTrkAux::fitTracksCm;
using AlephTrkAux::TrkBlock;
using AlephTrkAux::TrkAux;
using AlephTrkAux::pushTrk;
using AlephTrkAux::reserveTrk;

constexpr double M_K     = AlephMasses::kK;     // charged kaon
constexpr double M_PICH  = AlephMasses::kPiCh;  // charged pion
constexpr double M_D0    = AlephMasses::kD0;
constexpr double DM_NOMINAL = AlephMasses::kDeltaMDstarD0;

constexpr double E_BEAM = AlephUnits::kEBeam;  // for xE = E/E_BEAM
constexpr double PRE_MARGIN = 0.05;  // pre-fit mass window margin [GeV]

// Storage defaults, deliberately loose: the sample is re-cut offline.
constexpr double M_LO = 1.70, M_HI = 2.03;   // stored K pi mass window [GeV]
constexpr double CHI2_CUT = 25.;             // D0 vertex chi2 (ndf = 1)
constexpr double DPV_FID = 10.;              // |vtx-PV| storage fiducial [cm]
constexpr double DM_MAX = 0.20;              // stored dm ceiling [GeV]
// applied to the PERIGEE momentum, while the stored trk p is the at-vertex one
constexpr double P_MIN = 0.3;                // K/pi |p| floor [GeV]
constexpr double PS_MIN = 0.1;               // slow-pion |p| floor [GeV]

// Working-point flags, evaluated on the stored quantities.
constexpr double D0_LOOSE_DM = 0.060;        // |m(K pi) - m_D0| [GeV]
constexpr double D0_TIGHT_DM = 0.030;
constexpr double D0_TIGHT_DPVSIG = 3.0;      // |vtx-PV| significance
constexpr double D0_TIGHT_COSPOINT = 0.99;
constexpr double D0_TIGHT_COSSTAR = 0.8;     // |cos(theta*)| of the kaon
// D*: the dm window carries the tag
constexpr double DS_LOOSE_DM = 0.050, DS_LOOSE_DDM = 0.0030;
constexpr double DS_TIGHT_DM = 0.025, DS_TIGHT_DDM = 0.0015;
constexpr double DS_TIGHT_PS = 0.3;          // slow-pion |p| [GeV]
constexpr double DS_TIGHT_COSPOINT = 0.95;   // D0 vertex pointing
constexpr double TIGHT_PK = 1.0, TIGHT_PPI = 1.0;  // daughter |p| [GeV]
constexpr double TIGHT_CHI2 = 10.;           // D0 vertex chi2

inline TLorentzVector lorentz(const TVector3& p, double m) {
  return TLorentzVector(p, std::sqrt(p.Mag2() + m * m));
}

// Kaon helicity angle in the D0 rest frame, w.r.t. the D0 momentum direction in the lab.
inline double cosThetaStarK(const TVector3& pk, const TVector3& ppi) {
  TLorentzVector k = lorentz(pk, M_K);
  TLorentzVector d = k + lorentz(ppi, M_PICH);
  if (!(d.Vect().Mag() > 0.)) return -99.;
  TLorentzVector kr = k;
  kr.Boost(-d.BoostVector());
  double kk = kr.Vect().Mag();
  if (!(kk > 0.)) return -99.;
  return kr.Vect().Dot(d.Vect().Unit()) / kk;
}

// Kinematics block of the D0 and D* entries: p, px, py, pz, costheta and xE are
// those of the entry (a D* includes the slow pion); the rest are the D0's.
struct CandKin {
  RVec<float> m_kpi;           // K pi mass at the fitted vertex [GeV]
  RVec<float> p, px, py, pz, costheta, xE;
  RVec<float> chi2;            // vertex fit chi2, normalised (ndf = 1)
  RVec<float> vx, vy, vz;      // fitted vertex [cm]
  RVec<float> dpv, dpvSig;     // |vtx-PV| [cm] and its 3D significance
  RVec<float> cosPoint;        // cos(angle) between p(D0) and (vtx - PV)
  RVec<float> cosThetaStar;    // kaon helicity angle
};

// D0 -> K pi candidates, one entry per (pair, mass assignment).
struct D0Block {
  CandKin     kin;
  RVec<int>   loose, tight;    // labels, not cuts
  RVec<int>   nsec;            // legs in the secondary pool (0-2)
  TrkBlock    trkK, trkPi;
};

// D* -> D0 pi_slow candidates, one entry per (D0 candidate, third track).
struct DstarBlock {
  CandKin     kin;
  RVec<float> dm;
  RVec<int>   rs;              // 1 = right-sign slow pion (charge of the pi)
  RVec<int>   loose, tight;
  RVec<int>   d0idx;           // index of the parent entry in the internal D0 list
  RVec<int>   nsec;            // legs in the secondary pool (0-3)
  TrkBlock    trkK, trkPi, trkPis;
};

struct DstarCands {
  D0Block    d0;
  DstarBlock ds;
  int        nfits = 0;   // vertex fits actually performed
};

inline void reserveKin(CandKin& k, size_t n) {
  k.m_kpi.reserve(n); k.p.reserve(n); k.px.reserve(n); k.py.reserve(n);
  k.pz.reserve(n); k.costheta.reserve(n); k.xE.reserve(n); k.chi2.reserve(n);
  k.vx.reserve(n); k.vy.reserve(n); k.vz.reserve(n); k.dpv.reserve(n);
  k.dpvSig.reserve(n); k.cosPoint.reserve(n); k.cosThetaStar.reserve(n);
}

// One entry: p3 and energy of the entry itself, the rest of the D0 (pmag = |p3|).
inline void pushKin(CandKin& k, double m, const TVector3& p3, double pmag,
                    double energy, double chi2, const TVector3& x, double dpv,
                    float dpvSig, double cosPoint, double cosThetaStar) {
  k.m_kpi.push_back(m);
  k.p.push_back(pmag);
  k.px.push_back(p3.X());
  k.py.push_back(p3.Y());
  k.pz.push_back(p3.Z());
  k.costheta.push_back(p3.Z() / pmag);
  k.xE.push_back(energy / E_BEAM);
  k.chi2.push_back(chi2);
  k.vx.push_back(x.X());
  k.vy.push_back(x.Y());
  k.vz.push_back(x.Z());
  k.dpv.push_back(dpv);
  k.dpvSig.push_back(dpvSig);
  k.cosPoint.push_back(cosPoint);
  k.cosThetaStar.push_back(cosThetaStar);
}

// Reserve the output vectors for one entry per D0 candidate.
inline void reserveD0(D0Block& d, size_t n) {
  reserveKin(d.kin, n);
  d.loose.reserve(n); d.tight.reserve(n); d.nsec.reserve(n);
  reserveTrk(d.trkK, n); reserveTrk(d.trkPi, n);
}

inline void reserveDstar(DstarBlock& d, size_t n) {
  reserveKin(d.kin, n);
  d.dm.reserve(n); d.rs.reserve(n); d.loose.reserve(n);
  d.tight.reserve(n); d.d0idx.reserve(n); d.nsec.reserve(n);
  reserveTrk(d.trkK, n); reserveTrk(d.trkPi, n); reserveTrk(d.trkPis, n);
}

// The finder; every selection value is a constant defined above.
// The pool is the full baseline-selected list, primary and secondary alike.
// veto_orig = original-Tracks indices excluded from it.
inline DstarCands findDstar(
    const RVec<edm4hep::TrackState>& tracks,
    const RVec<int>& orig_idx,
    const RVec<int>& nvdet,
    const RVec<int>& nitc,
    const RVec<float>& chi2ndf,
    const RVec<int>& isprim,
    const RVec<int>& pool,
    const VertexingUtils::FCCAnalysesVertex& PV,
    const RVec<int>& veto_orig,
    double solenoidBz,
    double bsx, double bsy, double bsz) {

  DstarCands out;
  const int nTr = tracks.size();
  if (nTr < 2) return out;

  const TrkAux aux{&orig_idx, &nvdet, &nitc, &chi2ndf, &isprim, &pool, bsx, bsy, bsz};
  auto nv = [&](const RVec<int>& src, int k) {
    return (k >= 0 && k < (int)src.size()) ? src[k] : -1;
  };
  auto oidx = [&](int k) {
    return (k >= 0 && k < (int)orig_idx.size()) ? orig_idx[k] : -1;
  };
  // SEC = in the secondary set; "neither" tracks count as primary-like here
  auto isSec = [&](int k) { return k < (int)pool.size() && pool[k] == 1; };

  // momentum prefilter at the softest floor; canKPi marks the K/pi-capable
  // eK/ePi = hypothesis energies at the perigee, for the pre-fit mass window
  std::vector<int> good;
  good.reserve(nTr);
  std::vector<TVector3> pper(nTr);
  std::vector<double> eK(nTr, 0.), ePi(nTr, 0.);
  std::vector<TLorentzVector> pisLv(nTr);  // slow-pion four-momentum
  std::vector<char> canKPi(nTr, 0);
  const std::vector<char> vetoed = AlephTrkAux::memberMask(veto_orig);
  for (int i = 0; i < nTr; ++i) {
    const auto& t = tracks[i];
    if (AlephTrkAux::inMask(vetoed, oidx(i))) continue;
    pper[i] = perigeeMomentum(t, solenoidBz);
    double pm = pper[i].Mag();
    if (pm < PS_MIN) continue;
    canKPi[i] = (pm >= P_MIN) ? 1 : 0;
    eK[i] = std::sqrt(pper[i].Mag2() + M_K * M_K);
    ePi[i] = std::sqrt(pper[i].Mag2() + M_PICH * M_PICH);
    pisLv[i] = TLorentzVector(pper[i], ePi[i]);
    good.push_back(i);
  }
  if (good.size() < 2) return out;

  TVector3 pv(PV.vertex.position[0], PV.vertex.position[1], PV.vertex.position[2]);
  const double pre_lo = M_LO - PRE_MARGIN, pre_hi = M_HI + PRE_MARGIN;

  struct FitRes {
    bool ok;
    double chi2, dis;
    float dsig;
    TVector3 x, d, pa, pb;
  };
  RVec<edm4hep::TrackState> tr_pair(2);
  // i < j is guaranteed by the callers; pa/pb follow that order
  auto fitPair = [&](int i, int j) {
    FitRes r;
    r.ok = false;
    r.chi2 = 0.; r.dis = 0.; r.dsig = -1.;
    tr_pair[0] = tracks[i];
    tr_pair[1] = tracks[j];
    auto v = fitTracksCm(tr_pair, tracks, solenoidBz);
    ++out.nfits;
    if (v.updated_track_momentum_at_vertex.size() == 2) {
      double chi2 = v.vertex.chi2;
      bool pass = (chi2 == chi2) && !(chi2 >= CHI2_CUT);
      if (pass) {
        TVector3 x(v.vertex.position[0], v.vertex.position[1], v.vertex.position[2]);
        TVector3 dvec = x - pv;
        double dis = dvec.Mag();
        if (!(dis > DPV_FID)) {
          r.ok = true;
          r.chi2 = chi2; r.x = x; r.d = dvec; r.dis = dis;
          r.dsig = vertexDistSig(dvec, v.vertex.covMatrix, PV.vertex.covMatrix);
          r.pa = v.updated_track_momentum_at_vertex[0];
          r.pb = v.updated_track_momentum_at_vertex[1];
        }
      }
    }
    return r;
  };

  // one accepted (pair, mass assignment) before it is written out
  struct D0Cand {
    int iK, iPi;
    double m, chi2, dis, cosp, cstar, energy;
    float dsig;
    TVector3 p3, x, pK, pPi;
    int loose, tight;
  };

  // every opposite-charge pair inside the pre-fit window, both assignments
  auto makeD0 = [&]() {
    std::vector<D0Cand> res;
    for (size_t a = 0; a + 1 < good.size(); ++a) {
      const int i = good[a];
      if (!canKPi[i]) continue;
      for (size_t b = a + 1; b < good.size(); ++b) {
        const int j = good[b];
        if (!canKPi[j]) continue;
        if (tracks[i].omega * tracks[j].omega > 0) continue;  // opposite charge only

        // pre-fit window on the perigee momenta, both mass assignments
        const TVector3 p_pre = pper[i] + pper[j];
        const double p2_pre = p_pre.Mag2();
        const double e_ij = eK[i] + ePi[j], e_ji = ePi[i] + eK[j];
        double m_pre_ij = std::sqrt(std::max(0., e_ij * e_ij - p2_pre));
        double m_pre_ji = std::sqrt(std::max(0., e_ji * e_ji - p2_pre));
        bool ok_ij = (m_pre_ij >= pre_lo && m_pre_ij <= pre_hi);
        bool ok_ji = (m_pre_ji >= pre_lo && m_pre_ji <= pre_hi);
        if (!ok_ij && !ok_ji) continue;

        const FitRes fr = fitPair(i, j);
        if (!fr.ok) continue;

        for (int asg = 0; asg < 2; ++asg) {
          if (asg == 0 && !ok_ij) continue;
          if (asg == 1 && !ok_ji) continue;
          D0Cand c;
          c.iK  = (asg == 0) ? i : j;
          c.iPi = (asg == 0) ? j : i;
          c.pK  = (asg == 0) ? fr.pa : fr.pb;
          c.pPi = (asg == 0) ? fr.pb : fr.pa;

          c.m = AlephV0New::invMass(c.pK, M_K, c.pPi, M_PICH);
          if (c.m < M_LO || c.m > M_HI) continue;

          TLorentzVector d0lv = lorentz(c.pK, M_K) + lorentz(c.pPi, M_PICH);
          c.p3 = d0lv.Vect();
          double pmag = c.p3.Mag();
          if (!(pmag > 0.)) continue;
          c.energy = d0lv.E();
          c.chi2 = fr.chi2;
          c.x = fr.x;
          c.dis = fr.dis;
          c.dsig = fr.dsig;
          c.cosp = (fr.dis > 0.) ? c.p3.Dot(fr.d) / (pmag * fr.dis) : -99.;
          c.cstar = cosThetaStarK(c.pK, c.pPi);

          // labels evaluated on the STORED post-fit quantities
          const double dmass = std::abs(c.m - M_D0);
          c.loose = (dmass < D0_LOOSE_DM) ? 1 : 0;
          c.tight = (dmass < D0_TIGHT_DM && c.dsig > D0_TIGHT_DPVSIG &&
                     c.cosp > D0_TIGHT_COSPOINT &&
                     std::abs(c.cstar) < D0_TIGHT_COSSTAR &&
                     c.pK.Mag() > TIGHT_PK && c.pPi.Mag() > TIGHT_PPI &&
                     c.chi2 < TIGHT_CHI2) ? 1 : 0;
          res.push_back(c);
        }
      }
    }
    return res;
  };

  auto storeD0 = [&](const D0Cand& c) {
    pushKin(out.d0.kin, c.m, c.p3, c.p3.Mag(), c.energy, c.chi2, c.x, c.dis,
            c.dsig, c.cosp, c.cstar);
    out.d0.loose.push_back(c.loose);
    out.d0.tight.push_back(c.tight);
    out.d0.nsec.push_back((isSec(c.iK) ? 1 : 0) + (isSec(c.iPi) ? 1 : 0));
    pushTrk(out.d0.trkK, tracks, c.iK, c.pK, aux);
    pushTrk(out.d0.trkPi, tracks, c.iPi, c.pPi, aux);
  };

  // D* from a D0 candidate plus the slow-pion track k; no three-track fit
  auto pushDstar = [&](const D0Cand& c, int k, int d0idx) {
    // the D0 mass enters dm as the STORED (float) value
    const float mf = c.m;
    const double m = mf;
    TLorentzVector dslv = TLorentzVector(c.p3, c.energy) + pisLv[k];
    double dm = dslv.M() - m;
    if (!(dm > 0.)) return false;
    if (dm >= DM_MAX) return false;
    TVector3 p3 = dslv.Vect();
    double pmag = p3.Mag();
    if (!(pmag > 0.)) return false;

    const float chi2f = c.chi2;
    const int qPi = (tracks[c.iPi].omega > 0) ? 1 : -1;
    const int qs = (tracks[k].omega > 0) ? 1 : -1;

    pushKin(out.ds.kin, mf, p3, pmag, dslv.E(), chi2f, c.x, c.dis, c.dsig,
            c.cosp, c.cstar);
    out.ds.dm.push_back(dm);
    out.ds.rs.push_back((qs == qPi) ? 1 : 0);
    out.ds.d0idx.push_back(d0idx);
    out.ds.nsec.push_back((isSec(c.iK) ? 1 : 0) + (isSec(c.iPi) ? 1 : 0) +
                          (isSec(k) ? 1 : 0));
    pushTrk(out.ds.trkK, tracks, c.iK, c.pK, aux);
    pushTrk(out.ds.trkPi, tracks, c.iPi, c.pPi, aux);
    pushTrk(out.ds.trkPis, tracks, k, pper[k], aux);

    const double dmass = std::abs(m - M_D0);
    const double ddm = std::abs(dm - DM_NOMINAL);
    // primary-pattern veto: two primary-set D0 legs with a secondary slow pion
    const bool bad_pattern = (nv(isprim, c.iK) == 1 && nv(isprim, c.iPi) == 1 &&
                              nv(isprim, k) == 0);
    out.ds.loose.push_back((dmass < DS_LOOSE_DM && ddm < DS_LOOSE_DDM) ? 1 : 0);
    out.ds.tight.push_back((dmass < DS_TIGHT_DM && ddm < DS_TIGHT_DDM &&
                            c.pK.Mag() > TIGHT_PK && c.pPi.Mag() > TIGHT_PPI &&
                            chi2f < TIGHT_CHI2 && !bad_pattern &&
                            pper[k].Mag() > DS_TIGHT_PS &&
                            c.cosp > DS_TIGHT_COSPOINT) ? 1 : 0);
    return true;
  };

  // single all-track pass: every D0 candidate, then every third track
  auto cands = makeD0();
  reserveD0(out.d0, cands.size());
  reserveDstar(out.ds, cands.size());
  for (const auto& c : cands) storeD0(c);
  for (size_t c = 0; c < cands.size(); ++c)
    for (size_t s = 0; s < good.size(); ++s) {
      const int k = good[s];
      if (k == cands[c].iK || k == cands[c].iPi) continue;
      pushDstar(cands[c], k, (int)c);
    }
  return out;
}

// Staging class of each entry: 0 = primary set, 1 = secondary set, 2 = neither.
inline RVec<int> poolClass(const RVec<int>& orig_idx,
                           const RVec<int>& prim_orig,
                           const RVec<int>& sec_orig) {
  RVec<int> out;
  const std::vector<char> sec = AlephTrkAux::memberMask(sec_orig);
  const std::vector<char> prim = AlephTrkAux::memberMask(prim_orig);
  for (int o : orig_idx) {
    int v = 2;
    if (AlephTrkAux::inMask(sec, o)) v = 1;
    else if (AlephTrkAux::inMask(prim, o)) v = 0;
    out.push_back(v);
  }
  return out;
}

} // namespace AlephDstar
} // namespace FCCAnalyses

#endif
