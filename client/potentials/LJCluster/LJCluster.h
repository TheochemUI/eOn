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

// #include "../../system_unit.h" // unit converters
#include "client/Potential.h"
#include "client/potentials/LJ/LJ.h"
#include <limits>
namespace eonc {
/** Lennard Jones potential.*/
class LJCluster final : public Potential<LJCluster> {
  // NOTE(rg): Do we need LJCluster params?
private:
  LJ _lj;

public:
  LJCluster(const LJ::Params &p_a)
      : _lj{LJ(LJ::Params{.u0 = p_a.u0,
                          .cutoff_R = std::numeric_limits<double>::infinity(),
                          .psi = p_a.psi})} {}

  ~LJCluster();

  void forceImpl(const ForceInput &, ForceOut *) override final;
  void setParameters(double r0Recieved, double psiRecieved);
};
} // namespace eonc
