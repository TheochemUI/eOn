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
#include "../../Potential.h"

namespace eonc {
namespace pot {
struct LJParams {
  double u0{1.0};
  double cutoff{15.0};
  double psi{1.0};
  LJParams(const toml::table &tbl) {
    // TODO(rg)
    u0 = tbl["u0"].value_or(u0);
    cutoff = tbl["cutoff"].value_or(cutoff);
    psi = tbl["psi"].value_or(psi);
  }
};
} // namespace pot

/** Lennard Jones potential.*/
class LJ : public Potential<LJ> {
private:
  double u0;
  double cuttOffR;
  double psi;
  double cuttOffU;

public:
  LJ(pot::LJParams &ljp)
      : u0{ljp.u0},
        cuttOffR{ljp.cutoff},
        psi{ljp.psi} {}

  ~LJ() = default;

  void force(const ForceInput &params, ForceOut *efvdat) override;

  void setParameters(double r0Recieved, double u0Recieved, double psiRecieved);
};

} // namespace eonc
