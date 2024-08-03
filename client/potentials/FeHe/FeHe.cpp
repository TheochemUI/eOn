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
//-----------------------------------------------------------------------------------
#include "FeHe.h"
#include <assert.h>
#include <stdlib.h>
namespace eonc {

void FeHe::forceImpl(const ForceInput &fip, ForceOut *efvd) {
#ifdef EON_CHECKS
  eonc::pot::checkParams(fip);
  eonc::pot::zeroForceOut(fip.nAtoms, efvd);
#endif
  const long int N = fip.nAtoms;
  std::vector<double> RX(N), RY(N), RZ(N), FX(N), FY(N), FZ(N);
  std::vector<int> ISPEC(N);

  for (long int idx = 0; idx < N; ++idx) {
    RX[idx] = fip.pos[idx * 3 + 0];
    RY[idx] = fip.pos[idx * 3 + 1];
    RZ[idx] = fip.pos[idx * 3 + 2];
    assert(fip.atmnrs[idx] == 26 || fip.atmnrs[idx] == 2);
    ISPEC[idx] = (fip.atmnrs[idx] == 26) ? 0 : 1;
  }

  feforce_(&N, RX.data(), RY.data(), RZ.data(), ISPEC.data(), FX.data(),
           FY.data(), FZ.data(), &efvd->energy, &fip.box[0], &fip.box[4],
           &fip.box[8]);

  for (long int idx = 0; idx < N; ++idx) {
    efvd->F[idx * 3] = FX[idx];
    efvd->F[idx * 3 + 1] = FY[idx];
    efvd->F[idx * 3 + 2] = FZ[idx];
  }
}

} // namespace eonc
