#pragma once
#include "Matter.h"
#include <functional>

namespace eonc {

class StructComparer {
public:
  struct Params {
    ///< The distance criterion for comparing geometries
    double distanceDifference{0.1};
    ///< radius used in the local atomic structure analysis
    double neighborCutoff{3.3};
    bool checkRotation{false};
    bool indistinguishableAtoms{true};
    double energyDifference{0.01};
    bool removeTranslation{true};
  } scparams;

  // TODO(rg):: Indistinguishable conflicts with default Param
  bool compare(const Matter &m1, const Matter &m2,
               const bool indistinguishable = false);
  StructComparer(Params scp_a)
      : scparams{scp_a} {
    setupCompareFunc();
  }
  StructComparer(const toml::table &tbl) {
    fromTOML(tbl);
    setupCompareFunc();
  }

private:
  std::function<bool(const Matter &, const Matter &, double)> compareFunc;
  void fromTOML(const toml::table &tbl);
  void setupCompareFunc(); // Setup cutoff
};

} // namespace eonc