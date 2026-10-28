#ifndef ALEPH_ALEPHVARS_H
#define ALEPH_ALEPHVARS_H

/*
  ALEPH's own reconstruction results as exported by the converter, written as
  stage1 branches without any selection: the dE/dx hypotheses of every track
  (ALPHA QDEDX/QPADX, and QDEDXM in simulation), the best-match MC particle of
  every track (ALPHA KBESTM) and its parent, the YV0V V0 candidates, the
  beam spot and luminous region, the beam energy and the solenoid field.
  Lengths cm, momenta and masses GeV. kNotAvailable (-999)
  marks a value that is missing or must not be used.
*/

#include <cmath>
#include <stdexcept>
#include <string>

#include <ROOT/RVec.hxx>
#include "TVector3.h"

#include "edm4hep/MCParticleData.h"
#include "edm4hep/ParticleIDData.h"
#include "edm4hep/ReconstructedParticleData.h"
#include "edm4hep/VertexData.h"
#include "podio/ObjectID.h"

#include "aleph_units.h"

namespace FCCAnalyses {
namespace AlephVars {

using ROOT::VecOps::RVec;

constexpr float kNotAvailable = -999.f;
constexpr int kNotAvailableInt = -999;

// ---------------------------------------------------------------------------
// dE/dx hypotheses: dEdxWiresPID, dEdxPadsPID (13 parameters) and, in
// simulation, dEdxWiresQDEDXMPID (17), one ParticleID per Tracks entry
// (entry t belongs to Tracks[t]). Parameters: dEdx, nSamples, trackLength,
// then (expected_h, sigma_h) for h = electron, muon, pion, kaon, proton.
// ---------------------------------------------------------------------------
constexpr int kNHyp = 5;
constexpr int kNParDedx = 13;
constexpr int kNParQDEDXM = 17;

// Per-track arrays of one dE/dx collection. errorCode is the ParticleID type
// (ALPHA error code, 0 = measured) of every entry. A track without a usable
// measurement (type != 0, or dEdx <= 0) reads kNotAvailable in every other
// array; a hypothesis with ALPHA error code 5 (expected_h == 0) reads
// kNotAvailable in expected[h], sigma[h] and pull[h]. pull[h] = (value -
// expected[h]) / sigma[h], kNotAvailable if expected[h] is kNotAvailable or
// sigma[h] <= 0. With emptyIsMissing, an empty collection (QDEDXM in a
// simulated event without the ALPHA MC banks) reads kNotAvailable(Int) for
// every track.
struct TrackDedx {
  RVec<float> value, nSamples, trackLength;
  RVec<int> errorCode;
  RVec<float> expected[kNHyp], sigma[kNHyp], pull[kNHyp];
};

inline TrackDedx trackDedx(const RVec<edm4hep::ParticleIDData>& pid,
                           const RVec<float>& params, size_t nTracks,
                           int nPar, const std::string& name,
                           bool emptyIsMissing = false) {
  // index-parallel layout: a different size would attribute every value to the
  // wrong track, so the input is refused
  if (pid.size() != nTracks && !(emptyIsMissing && pid.empty()))
    throw std::runtime_error("AlephVars::trackDedx: " + name + " has " +
                             std::to_string(pid.size()) + " entries for " +
                             std::to_string(nTracks) +
                             " Tracks; it must hold exactly one ParticleID per track");
  TrackDedx out;
  out.value = RVec<float>(nTracks, kNotAvailable);
  out.nSamples = RVec<float>(nTracks, kNotAvailable);
  out.trackLength = RVec<float>(nTracks, kNotAvailable);
  out.errorCode = RVec<int>(nTracks, kNotAvailableInt);
  for (int h = 0; h < kNHyp; ++h) {
    out.expected[h] = RVec<float>(nTracks, kNotAvailable);
    out.sigma[h] = RVec<float>(nTracks, kNotAvailable);
    out.pull[h] = RVec<float>(nTracks, kNotAvailable);
  }
  for (size_t t = 0; t < pid.size(); ++t) {
    const auto& e = pid[t];
    if (e.parameters_end < e.parameters_begin ||
        e.parameters_end - e.parameters_begin != (unsigned)nPar ||
        e.parameters_end > params.size())
      throw std::runtime_error("AlephVars::trackDedx: " + name + " entry " +
                               std::to_string(t) + " does not hold " +
                               std::to_string(nPar) + " parameters");
    const float* p = params.data() + e.parameters_begin;
    out.errorCode[t] = e.type;
    if (e.type != 0 || !(p[0] > 0.f)) continue;
    out.value[t] = p[0];
    out.nSamples[t] = p[1];
    out.trackLength[t] = p[2];
    for (int h = 0; h < kNHyp; ++h) {
      if (p[3 + 2 * h] == 0.f) continue;  // error code 5 for this hypothesis
      out.expected[h][t] = p[3 + 2 * h];
      out.sigma[h][t] = p[4 + 2 * h];
      if (p[3 + 2 * h] != kNotAvailable && p[4 + 2 * h] > 0.f)
        out.pull[h][t] = (p[0] - p[3 + 2 * h]) / p[4 + 2 * h];
    }
  }
  return out;
}

// Truth of every track from ALPHA's best match (KBESTM), the
// dEdxWiresQDEDXMMCLink relation (at most one link per track): the PDG code of
// the best-match MC particle and of its immediate parent in the MCParticles
// parents relation (0 = no best match or no parent), the MCParticles index of
// both (-1 = none), and the momentum and production vertex (transverse radius
// and z) of the best-match particle (kNotAvailable = no best match). From
// trackMCLink (track <-> MC particle, weight = fraction of shared track
// coordinates): the number of links of the track and the weight of its link to
// the best-match particle (the first such link; -1 = not linked there,
// kNotAvailable = no best match). noQDEDXM (empty QDEDXM collections: a
// simulated event without the ALPHA MC banks) gives kNotAvailable(Int) for
// every track; the link must then be empty too.
struct TrackTruth {
  RVec<int> pdg, motherPdg, idx, motherIdx;
  RVec<float> p, vtxR, vtxZ;
  RVec<int> nLinks;
  RVec<float> sharedFrac;
};

inline TrackTruth trackBestMatchTruth(size_t nTracks, bool noQDEDXM,
                                      const RVec<podio::ObjectID>& from,
                                      const RVec<podio::ObjectID>& to,
                                      const RVec<edm4hep::MCParticleData>& mc,
                                      const RVec<podio::ObjectID>& mcParents,
                                      const RVec<podio::ObjectID>& linkFrom,
                                      const RVec<podio::ObjectID>& linkTo,
                                      const RVec<float>& linkWeight) {
  if (noQDEDXM && !from.empty())
    throw std::runtime_error("AlephVars::trackBestMatchTruth: dEdxWiresQDEDXMMCLink has " +
                             std::to_string(from.size()) +
                             " entries but dEdxWiresQDEDXMPID is empty");
  const RVec<float> na(nTracks, kNotAvailable);
  const RVec<int> naInt(nTracks, kNotAvailableInt);
  if (noQDEDXM) return {naInt, naInt, naInt, naInt, na, na, na, naInt, na};
  if (from.size() != to.size())
    throw std::runtime_error("AlephVars::trackBestMatchTruth: link size mismatch");
  TrackTruth out{RVec<int>(nTracks, 0), RVec<int>(nTracks, 0), RVec<int>(nTracks, -1),
                 RVec<int>(nTracks, -1), na, na, na, RVec<int>(nTracks, 0), na};
  RVec<char> seen(nTracks, 0);
  for (size_t i = 0; i < from.size(); ++i) {
    const int t = from[i].index, m = to[i].index;
    if (t < 0 || (size_t)t >= nTracks || m < 0 || (size_t)m >= mc.size())
      throw std::runtime_error(
          "AlephVars::trackBestMatchTruth: dEdxWiresQDEDXMMCLink entry out of range");
    if (seen[t])
      throw std::runtime_error(
          "AlephVars::trackBestMatchTruth: two dEdxWiresQDEDXMMCLink entries for one track");
    seen[t] = 1;
    const auto& p = mc[m];
    out.pdg[t] = p.PDG;
    out.idx[t] = m;
    out.sharedFrac[t] = -1.f;
    out.p[t] = std::sqrt((double)p.momentum.x * p.momentum.x + (double)p.momentum.y * p.momentum.y +
                         (double)p.momentum.z * p.momentum.z);
    out.vtxR[t] = std::hypot((double)p.vertex.x, (double)p.vertex.y);
    out.vtxZ[t] = p.vertex.z;
    if (p.parents_end <= p.parents_begin) continue;
    const int q = p.parents_end <= mcParents.size() ? mcParents[p.parents_begin].index : -1;
    if (q < 0 || (size_t)q >= mc.size())
      throw std::runtime_error("AlephVars::trackBestMatchTruth: MCParticles parent out of range");
    out.motherPdg[t] = mc[q].PDG;
    out.motherIdx[t] = q;
  }
  if (linkFrom.size() != linkTo.size() || linkWeight.size() != linkFrom.size())
    throw std::runtime_error("AlephVars::trackBestMatchTruth: trackMCLink size mismatch");
  for (size_t i = 0; i < linkFrom.size(); ++i) {
    const int t = linkFrom[i].index, m = linkTo[i].index;
    if (t < 0 || (size_t)t >= nTracks || m < 0 || (size_t)m >= mc.size())
      throw std::runtime_error("AlephVars::trackBestMatchTruth: trackMCLink entry out of range");
    ++out.nLinks[t];
    if (m == out.idx[t] && out.sharedFrac[t] == -1.f) out.sharedFrac[t] = linkWeight[i];
  }
  return out;
}

// ---------------------------------------------------------------------------
// ALEPH V0 candidates (YV0V): V0Vertices[i] (position, covariance, chi2,
// relation to the daughters' RecoParticles) and V0VerticesRefitted[i] (30
// YV0V words as parameters, relation to the refitted daughters P1 = positive,
// P2 = negative in RefittedParticles).
// ---------------------------------------------------------------------------
constexpr int kNParV0 = 30;

struct AlephV0s {
  RVec<float> vx, vy, vz, cov[6], chi2, helixDistance, alpha, qt, covAlphaAlpha,
      covAlphaQt, covQtQt, px, py, pz, mKs, mLambda, mLambdaBar, mGamma;
  RVec<int> fitCode, trkPos, trkNeg;
};

inline float pairMass(const TVector3& a, double ma, const TVector3& b, double mb) {
  const double e = std::sqrt(a.Mag2() + ma * ma) + std::sqrt(b.Mag2() + mb * mb);
  const double m2 = e * e - (a + b).Mag2();
  return m2 > 0. ? std::sqrt(m2) : 0.f;
}

// requireWords = false accepts V0VerticesRefitted without the YV0V words
// (inputs converted before they were exported): the words read kNotAvailable.
// The daughter tracks are found through V0Vertices.particles -> RecoParticles
// -> Tracks; the converter links a daughter only when it has an energy-flow
// object, in RecoParticles order, so the charge of that object tells the
// positive from the negative daughter. -1 = no linked track.
inline AlephV0s alephV0s(const RVec<edm4hep::VertexData>& v0,
                         const RVec<podio::ObjectID>& v0Particles,
                         const RVec<edm4hep::VertexData>& v0r,
                         const RVec<float>& v0rParams,
                         const RVec<podio::ObjectID>& v0rParticles,
                         const RVec<edm4hep::ReconstructedParticleData>& refitted,
                         const RVec<edm4hep::ReconstructedParticleData>& rps,
                         const RVec<podio::ObjectID>& rpTracks, size_t nTracks,
                         bool requireWords) {
  if (v0r.size() != v0.size())
    throw std::runtime_error("AlephVars::alephV0s: V0VerticesRefitted and V0Vertices differ in size");
  AlephV0s out;
  const auto bad = [](const std::string& what) {
    return std::runtime_error("AlephVars::alephV0s: " + what);
  };
  for (size_t j = 0; j < v0.size(); ++j) {
    const auto& v = v0[j];
    out.vx.push_back(v.position.x);
    out.vy.push_back(v.position.y);
    out.vz.push_back(v.position.z);
    for (int c = 0; c < 6; ++c) out.cov[c].push_back(v.covMatrix.values[c]);
    out.chi2.push_back(v.chi2);

    // YV0V words
    const auto& r = v0r[j];
    const unsigned np = r.parameters_end >= r.parameters_begin ? r.parameters_end - r.parameters_begin : 0;
    if (np == (unsigned)kNParV0 && r.parameters_end <= v0rParams.size()) {
      const float* w = v0rParams.data() + r.parameters_begin;
      out.alpha.push_back(w[0]);
      out.qt.push_back(w[1]);
      out.covAlphaAlpha.push_back(w[2]);
      out.covAlphaQt.push_back(w[3]);
      out.covQtQt.push_back(w[4]);
      out.fitCode.push_back((int)std::lround(w[5]));
      out.helixDistance.push_back(w[6]);
    } else if (np == 0 && !requireWords) {
      for (auto* a : {&out.alpha, &out.qt, &out.covAlphaAlpha, &out.covAlphaQt,
                      &out.covQtQt, &out.helixDistance})
        a->push_back(kNotAvailable);
      out.fitCode.push_back(kNotAvailableInt);
    } else {
      throw bad("V0VerticesRefitted entry " + std::to_string(j) + " holds " +
                std::to_string(np) + " parameters instead of the " +
                std::to_string(kNParV0) + " YV0V words" +
                (np == 0 ? " (input converted before they were exported)" : ""));
    }

    // refitted daughters: momentum sum and the pair masses of the hypotheses
    if (r.particles_end - r.particles_begin != 2 || r.particles_end > v0rParticles.size())
      throw bad("V0VerticesRefitted entry " + std::to_string(j) + " does not link two refitted daughters");
    TVector3 pd[2];
    float qd[2];
    for (int k = 0; k < 2; ++k) {
      const int i = v0rParticles[r.particles_begin + k].index;
      if (i < 0 || (size_t)i >= refitted.size()) throw bad("refitted daughter out of range");
      const auto& m = refitted[i].momentum;
      pd[k].SetXYZ(m.x, m.y, m.z);
      qd[k] = refitted[i].charge;
    }
    if (!(qd[0] > 0.f && qd[1] < 0.f)) throw bad("refitted daughters are not (positive, negative)");
    const TVector3 sum = pd[0] + pd[1];
    out.px.push_back(sum.X());
    out.py.push_back(sum.Y());
    out.pz.push_back(sum.Z());
    using namespace AlephMasses;
    out.mKs.push_back(pairMass(pd[0], kPiCh, pd[1], kPiCh));
    out.mLambda.push_back(pairMass(pd[0], kProton, pd[1], kPiCh));
    out.mLambdaBar.push_back(pairMass(pd[0], kPiCh, pd[1], kProton));
    out.mGamma.push_back(pairMass(pd[0], kElectron, pd[1], kElectron));

    // daughter tracks
    int trk[2] = {-1, -1};  // positive, negative
    if (v.particles_end > v0Particles.size()) throw bad("V0Vertices particle relation out of range");
    for (unsigned k = v.particles_begin; k < v.particles_end; ++k) {
      const int ip = v0Particles[k].index;
      if (ip < 0 || (size_t)ip >= rps.size()) throw bad("daughter RecoParticle out of range");
      const auto& rp = rps[ip];
      if (rp.tracks_end <= rp.tracks_begin) continue;
      const int t = rp.tracks_begin < rpTracks.size() ? rpTracks[rp.tracks_begin].index : -1;
      if (t < 0 || (size_t)t >= nTracks) throw bad("daughter track out of range");
      if (rp.charge == 0.f) throw bad("daughter RecoParticle without charge");
      const int s = rp.charge > 0.f ? 0 : 1;
      if (trk[s] >= 0) throw bad("two daughters of the same charge");
      trk[s] = t;
    }
    out.trkPos.push_back(trk[0]);
    out.trkNeg.push_back(trk[1]);
  }
  return out;
}

// ---------------------------------------------------------------------------
// Beam spot (BeamSpot) and luminous region (LuminousRegion), at most one entry
// each per event; the uncertainties and sizes are the square roots of the
// stored covariance diagonal (the off-diagonal terms are zero).
// ---------------------------------------------------------------------------
struct AlephBeamSpot {
  float x = kNotAvailable, y = kNotAvailable, z = kNotAvailable;
  float errx = kNotAvailable, erry = kNotAvailable, errz = kNotAvailable;
  int algo = kNotAvailableInt;
  float sizex = kNotAvailable, sizey = kNotAvailable, sizez = kNotAvailable;
};

inline float sqrtOrNA(float v) { return v >= 0.f ? std::sqrt(v) : kNotAvailable; }

inline AlephBeamSpot alephBeamSpot(const RVec<edm4hep::VertexData>& bs,
                                   const RVec<edm4hep::VertexData>& lr) {
  if (bs.size() > 1 || lr.size() > 1)
    throw std::runtime_error("AlephVars::alephBeamSpot: more than one BeamSpot or LuminousRegion entry in an event");
  AlephBeamSpot out;
  if (bs.size() == 1) {
    const auto& b = bs[0];
    out.x = b.position.x;
    out.y = b.position.y;
    out.z = b.position.z;
    out.errx = sqrtOrNA(b.covMatrix.values[0]);
    out.erry = sqrtOrNA(b.covMatrix.values[2]);
    out.errz = sqrtOrNA(b.covMatrix.values[5]);
    out.algo = b.algorithmType;
  }
  if (lr.size() == 1) {
    out.sizex = sqrtOrNA(lr[0].covMatrix.values[0]);
    out.sizey = sqrtOrNA(lr[0].covMatrix.values[2]);
    out.sizez = sqrtOrNA(lr[0].covMatrix.values[5]);
  }
  return out;
}

// ---------------------------------------------------------------------------
// RunInformation (one float list per event), filled by the converter in this
// order, each group only if its bank exists: RLEP LEP beam energy [MeV], MBFL
// solenoid field [kG], RALE magnet polarity and three currents, MBLC LCAL and
// SICAL luminosity. The size tells the groups apart: 8 or 4 (data / simulation
// without RALE) hold the field, 7 or 3 come from a converter before the field
// was exported (requireField = false accepts them: the field reads
// kNotAvailable). Any other size stops stage1, and so does a first value that
// is not a beam energy in MeV: a list without RLEP can have these sizes too.
// ---------------------------------------------------------------------------
struct AlephRunInfo {
  float beamEnergy = kNotAvailable, bz = kNotAvailable;  // GeV, T
};

// The first RunInformation value must be a LEP beam energy in MeV; a field
// in kG or a polarity is far below it.
constexpr float kMinBeamEnergyMeV = 1000.f;

inline AlephRunInfo alephRunInfo(const RVec<float>& ri, bool requireField) {
  AlephRunInfo out;
  const size_t n = ri.size();
  const bool withField = n == 8 || n == 4;
  if (!withField && !((n == 7 || n == 3) && !requireField))
    throw std::runtime_error("AlephVars::alephRunInfo: RunInformation holds " + std::to_string(n) +
                             " values, not the 8 or 4 of a converter that exports the magnetic field" +
                             (n == 7 || n == 3 ? " (older converter: pass --oldConverter)" : ""));
  if (!(ri[0] > kMinBeamEnergyMeV))
    throw std::runtime_error("AlephVars::alephRunInfo: the first of the " + std::to_string(n) +
                             " RunInformation values, " + std::to_string(ri[0]) +
                             ", is not a LEP beam energy in MeV (RLEP bank missing?)");
  out.beamEnergy = ri[0] / 1000.f;
  if (withField) out.bz = ri[1] / 10.f;
  return out;
}

}  // namespace AlephVars
}  // namespace FCCAnalyses

#endif  // ALEPH_ALEPHVARS_H
