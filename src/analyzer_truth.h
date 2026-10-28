#ifndef ALEPHTRUTH_H
#define ALEPHTRUTH_H

/*
  The secondary-track index recovery (secondaryToOriginalTrack) and the
  event-order V0 candidate getters. Nothing here reads MC truth: everything
  runs on data too. Positions in cm.
*/

#include <stdexcept>

#include <ROOT/RVec.hxx>
#include "TVector3.h"

#include "edm4hep/TrackState.h"

#include "FCCAnalyses/VertexingUtils.h"

#include "analyzer_trkaux.h"

namespace FCCAnalyses {
namespace AlephTruth {

using ROOT::VecOps::RVec;

// V0 candidates hold reco_ind into SecondaryTracks_looseBS: for each secondary
// (flipped param space) find its position in the flipped selected-baseline list
// and return its index in the original Tracks collection.
inline RVec<int> secondaryToOriginalTrack(
    const RVec<edm4hep::TrackState>& secondaries,
    const RVec<edm4hep::TrackState>& selected_flipped,
    const RVec<int>& selected_orig_idx) {
  RVec<int> out;
  out.reserve(secondaries.size());
  for (const auto& s : secondaries) {
    int found = -1;
    for (size_t k = 0; k < selected_flipped.size(); ++k) {
      if (VertexingUtils::compare_Tracks(s, selected_flipped[k])) { found = selected_orig_idx[k]; break; }
    }
    if (found < 0)
      throw std::runtime_error("AlephTruth: secondary track not found in selected baseline");
    out.push_back(found);
  }
  return out;
}

// Event-order candidate kinematics (independent of jet assignment). Same
// quantities as the VertexingUtils get_d3d_SV / get_chi2_SV / get_pMag_SV /
// get_x_SV getters, but RVec<float> rather than RVec<double>; candChi2 is the
// normalised chi2 (get_norm_chi2_SV), equal to get_chi2_SV only for a
// two-track vertex (ndf = 1).
// candDxyz returns kUndef when the reference vertex has fewer than
// kPVMinTracks tracks, i.e. is the default vertex at the origin.
inline RVec<float> candDxyz(const VertexingUtils::FCCAnalysesV0& v0s,
                            const VertexingUtils::FCCAnalysesVertex& PV) {
  RVec<float> out;
  if (PV.ntracks < AlephTrkAux::kPVMinTracks)
    return RVec<float>(v0s.vtx.size(), AlephTrkAux::kUndef);
  TVector3 pv(PV.vertex.position[0], PV.vertex.position[1], PV.vertex.position[2]);
  for (const auto& v : v0s.vtx) {
    TVector3 x(v.vertex.position[0], v.vertex.position[1], v.vertex.position[2]);
    out.push_back((x - pv).Mag());
  }
  return out;
}

// Fitted-vertex position component, axis 0/1/2 = x/y/z, in cm.
inline RVec<float> candVtxPos(const VertexingUtils::FCCAnalysesV0& v0s, int axis) {
  RVec<float> out;
  for (const auto& v : v0s.vtx) out.push_back(v.vertex.position[axis]);
  return out;
}

inline RVec<float> candChi2(const VertexingUtils::FCCAnalysesV0& v0s) {
  RVec<float> out;
  for (const auto& v : v0s.vtx) out.push_back(v.vertex.chi2);
  return out;
}

inline RVec<float> candP(const VertexingUtils::FCCAnalysesV0& v0s) {
  RVec<float> out;
  for (const auto& v : v0s.vtx)
    out.push_back(AlephTrkAux::candMomentum(v).Mag());
  return out;
}

// Component comp (0/1/2 = x/y/z) of the SAME summed vertex momentum whose
// magnitude candP returns, so sqrt(px^2+py^2+pz^2) reproduces candP exactly.
inline RVec<float> candPcomp(const VertexingUtils::FCCAnalysesV0& v0s, int comp) {
  RVec<float> out;
  for (const auto& v : v0s.vtx)
    out.push_back(AlephTrkAux::candMomentum(v)[comp]);
  return out;
}

inline RVec<float> candCosPointing(const VertexingUtils::FCCAnalysesV0& v0s,
                                   const VertexingUtils::FCCAnalysesVertex& PV) {
  RVec<float> out;
  TVector3 pv(PV.vertex.position[0], PV.vertex.position[1], PV.vertex.position[2]);
  for (const auto& v : v0s.vtx) {
    TVector3 x(v.vertex.position[0], v.vertex.position[1], v.vertex.position[2]);
    TVector3 p = AlephTrkAux::candMomentum(v);
    TVector3 d = x - pv;
    out.push_back((d.Mag() > 0 && p.Mag() > 0) ? d.Dot(p) / (d.Mag() * p.Mag()) : -2.);
  }
  return out;
}

} // namespace AlephTruth
} // namespace FCCAnalyses

#endif
