#ifndef ALEPH_UNITS_H
#define ALEPH_UNITS_H
// Unit conventions shared by the ALEPH analyzers: lengths cm, momenta GeV, field T.
// Dependency-free, so that standalone builds can include it as well.
namespace FCCAnalyses {
namespace AlephUnits {
// pT [GeV] = kPtPerTeslaCm * Bz [T] / |omega [1/cm]|  (0.29979 GeV/(T m), in cm)
constexpr double kPtPerTeslaCm = 0.0029979;
}  // namespace AlephUnits
}  // namespace FCCAnalyses
#endif
