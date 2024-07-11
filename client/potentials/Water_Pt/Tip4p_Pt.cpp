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
void Tip4p_Pt::force(long N, const double *R, const int *atomicNrs, double *F,
                     double *U, double *variance, const double *box) {
  variance = nullptr;
  double diagbox[3];
  diagbox[0] = box[0];
  diagbox[1] = box[4];
  diagbox[2] = box[8];
  int i = 0;
  while (atomicNrs[i] == 1)
    i += 2;
  computeHH_O_Pt_(i / 2, N - i * 3 / 2, R, F, *U, diagbox, 0);
}

} // namespace eonc
