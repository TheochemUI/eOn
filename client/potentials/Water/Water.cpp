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
#include "Water.hpp"

void Tip4p::force(long N, const double *R, const int *atomicNrs, double *F,
                  double *U, const double *box) {
  double diagbox[3];
  diagbox[0] = box[0];
  diagbox[1] = box[4];
  diagbox[2] = box[8];
  computeHH_O_(N, R, F, *U, diagbox, 0);
}

void SpceCcl::force(long N, const double *R, const int *atomicNrs, double *F,
                    double *U, double *variance, const double *box) {
  variance = nullptr;
  double diagbox[3];
  diagbox[0] = box[0];
  diagbox[1] = box[4];
  diagbox[2] = box[8];
  computeHH_O_(N, R, F, *U, diagbox, 0);
}
