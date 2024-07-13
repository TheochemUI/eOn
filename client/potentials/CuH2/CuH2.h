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

#pragma once
#include "client/Potential.h"

// natms(2), ndim, U(1), R(ndim), F(ndim), box(3)
extern "C" {
void c_force_eam(int natms[2], int ndim, double box[3], double *R, double *F,
                 double U[1]);
}
namespace eonc {
class CuH2 final : public Potential<CuH2> {

public:
  // NOTE(miha) you can recenter the coordinates before and after the calls
  void forceImpl(const ForceInput &, ForceOut *) override final;
};

} // namespace eonc
