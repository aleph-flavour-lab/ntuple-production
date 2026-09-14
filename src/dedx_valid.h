#ifndef ALEPH_DEDX_VALID_H
#define ALEPH_DEDX_VALID_H

#include <cmath>

namespace FCCAnalyses {
namespace AlephDedx {

// Single source of the dE/dx validity rule: a failed leg copies the track's omega into its
// value, so valid = finite positive value/error and |value| != |omega| (sign-convention free).
inline bool dEdxValid(float value, float error, float trackOmega) {
  return std::isfinite(value) && value > 0.f &&
         std::abs(value) != std::abs(trackOmega) && std::isfinite(error) &&
         error > 0.f;
}

} // namespace AlephDedx
} // namespace FCCAnalyses

#endif
