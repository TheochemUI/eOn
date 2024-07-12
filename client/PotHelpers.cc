#include "PotHelpers.hpp"

namespace eonc::pot {

// Typically this is done by the caller, however it is here as a sanity check
void zeroForceOut(const size_t &nAtoms, ForceOut *efvd) {
  efvd->energy = 0;
  efvd->variance = 0;
  for (size_t idx{0}; idx < nAtoms; idx++) {
    efvd->F[3 * idx] = 0;
    efvd->F[3 * idx + 1] = 0;
    efvd->F[3 * idx + 2] = 0;
  }
};

} // namespace eonc::pot
