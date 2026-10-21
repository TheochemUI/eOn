#pragma once

#include "rgpot/units.hpp"

namespace cpmd_stress_oracle {

inline constexpr double kGradientHartree = 0.5;

inline constexpr double nativeStress(int component) {
  return 0.01 * (static_cast<double>(component) + 1.0);
}

inline double stressEvPerAngstrom3(int component) {
  return nativeStress(component) *
         rgpot::units::HARTREE_PER_BOHR3_TO_EV_PER_ANGSTROM3;
}

inline double sessionHartree(double cell_zz) {
  return (0.75 + 0.001 * cell_zz) / rgpot::units::HARTREE_TO_EV;
}

} // namespace cpmd_stress_oracle
