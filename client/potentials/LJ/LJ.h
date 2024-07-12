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

/** Lennard Jones potential.*/
class LJ final : public Potential<LJ> {
public:
  struct Params {
    double u0{1.0};
    double cutoff_R{15.0};
    double psi{1.0};
    Params(const toml::table &tbl) {
      u0 = tbl["u0"].value_or(u0);
      cutoff_R = tbl["cutoff"].value_or(cutoff_R);
      psi = tbl["psi"].value_or(psi);
    }
  };

private:
  double u0;
  double cutoff_R;
  double psi;
  double cutoff_U{0};

  double calc_cutoffU(const Params &p);

public:
  LJ(const Params &ljp)
      : u0{ljp.u0},
        cutoff_R{ljp.cutoff_R},
        psi{ljp.psi},
        cutoff_U{calc_cutoffU(ljp)} {}

  void forceImpl(const ForceInput &params, ForceOut *efvdat) override;

  void setParameters(const Params &ljp);
};

} // namespace eonc
