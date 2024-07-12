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

namespace eonc {

/** Lennard Jones potential.*/
class LJ final : public Potential<LJ> {
public:
  struct Params {
    double u0{1.0};
    double cutoff_R{15.0};
    double psi{1.0};
    Params(const toml::table &tbl) {
      u0 = tbl["Potential"]["LJ"]["u0"].value_or(u0);
      cutoff_R = tbl["Potential"]["LJ"]["cutoff"].value_or(cutoff_R);
      psi = tbl["Potential"]["LJ"]["psi"].value_or(psi);
    }
  };

private:
  double u0;
  double cutoff_R;
  double psi;
  double cutoff_U{0};

  double calc_cutoffU(const LJ::Params &);

public:
  LJ(const LJ::Params &p_a)
      : u0{p_a.u0},
        cutoff_R{p_a.cutoff_R},
        psi{p_a.psi},
        cutoff_U{calc_cutoffU(p_a)} {}

  void forceImpl(const ForceInput &, ForceOut *) override final;

  void setParameters(const LJ::Params &);
};

} // namespace eonc
