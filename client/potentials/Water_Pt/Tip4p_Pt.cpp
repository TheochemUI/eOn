/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
**
** Copyright (c) 2010--present, eOn Development Team
** All rights reserved.
**
** Repo:
** https://github.com/TheochemUI/eOn
*/
#include "Tip4p_Pt.hpp"
namespace eonc {
void Tip4p_Pt::forceImpl(const ForceInput &fip, ForceOut *efvd) {
#ifdef EON_CHECKS
  eonc::pot::checkParams(fip);
  eonc::pot::zeroForceOut(fip.nAtoms, efvd);
#endif
  const long int N = fip.nAtoms;
  std::array<double, 3> diagbox{fip.box[0], fip.box[4], fip.box[8]};
  int idx{0};
  while (fip.atmnrs[idx] == 1) {
    idx += 2;
  }
  computeHH_O_Pt_(idx / 2, N - idx * 3 / 2, fip.pos, efvd->F, efvd->energy,
                  diagbox.data(), 0);
}

} // namespace eonc
