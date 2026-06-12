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
#include "FeHe.h"
#include <cassert>
#include <vector>

void FeHe::force(long N, const double *R, const int *atomicNrs, double *F,
                 double *U, double *variance, const double *box) {
  variance = nullptr;
  assert(N > 1);

  std::vector<double> RX(N), RY(N), RZ(N);
  std::vector<double> FX(N), FY(N), FZ(N);
  std::vector<int> ISPEC(N);

  for (long i = 0; i < N; i++) {
    RX[i] = R[i * 3 + 0];
    RY[i] = R[i * 3 + 1];
    RZ[i] = R[i * 3 + 2];
    assert(atomicNrs[i] == 26 || atomicNrs[i] == 2);
    ISPEC[i] = (atomicNrs[i] == 26) ? 0 : 1;
  }

  m_feforce(&N, RX.data(), RY.data(), RZ.data(), ISPEC.data(), FX.data(),
            FY.data(), FZ.data(), U, &box[0], &box[4], &box[8]);

  for (long i = 0; i < N; i++) {
    F[i * 3] = FX[i];
    F[i * 3 + 1] = FY[i];
    F[i * 3 + 2] = FZ[i];
  }
}
