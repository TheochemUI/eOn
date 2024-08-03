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
#include "../../Potential.h"
namespace eonc {
class ExtPot : public Potential<ExtPot> {

public:
  ExtPot(std::string extPotPath)
      : eon_extpot_path{extPotPath.c_str()} {};
  void forceImpl(const ForceInput &, ForceOut *) override final;

private:
  void passToSystem(long N, const double *R, const size_t *atomicNrs,
                    const double *box);
  void receiveFromSystem(long N, double *F, double *U);
  std::string eon_extpot_path;
};

} // namespace eonc
