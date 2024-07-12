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
#include <cmath>
namespace eonc {
class Morse final : public Potential<Morse> {
public:
  struct Params {
    // De in eV, a in Angstroms, re in Angstroms, cutoff in Angstroms
    double De{0.7102};
    double a{1.6047};
    double re{2.8970};
    double cutoff{9.5};
    Params(const toml::table &tbl) {
      De = tbl["Potential"]["De"].value_or(De);
      a = tbl["Potential"]["a"].value_or(a);
      re = tbl["Potential"]["re"].value_or(re);
      cutoff = tbl["Potential"]["cutoff"].value_or(cutoff);
    }
  };

private:
  void morse(double r, double &energy, double &force);
  double De_;
  double a_;
  double re_;
  double cutoff_;
  double energyCutoff_;

public:
  Morse(const Params &mpar)
      : De_{mpar.De},
        a_{mpar.a},
        re_{mpar.re},
        cutoff_{mpar.cutoff} {}
  void forceImpl(const ForceInput &params, ForceOut *efvdat) override final;
  void setParameters(const Params &mpar);
};

} // namespace eonc
