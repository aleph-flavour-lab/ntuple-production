#ifndef ALEPH_UNITS_H
#define ALEPH_UNITS_H
// Unit conventions shared by the ALEPH analyzers: lengths cm, momenta GeV, field T.
// Dependency-free, so that standalone builds can include it as well.
namespace FCCAnalyses {
namespace AlephUnits {
// pT [GeV] = kPtPerTeslaCm * Bz [T] / |omega [1/cm]|  (0.29979 GeV/(T m), in cm)
constexpr double kPtPerTeslaCm = 0.0029979;
}  // namespace AlephUnits

// PDG 2024 central values [GeV]: the single source for every analyzer mass.
// The charged-pion constant cannot be called M_PI, that name is a <cmath> macro.
namespace AlephMasses {
constexpr double kPiCh   = 0.13957039;
constexpr double kProton = 0.93827208;
constexpr double kKs     = 0.497611;
constexpr double kLambda = 1.115683;
}  // namespace AlephMasses
}  // namespace FCCAnalyses
#endif
