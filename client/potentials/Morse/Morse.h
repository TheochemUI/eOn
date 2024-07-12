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
    Params(const toml::node_view<const toml::node> &tbl) {
      De = tbl["Potential"]["Morse"]["De"].value_or(De);
      a = tbl["Potential"]["Morse"]["a"].value_or(a);
      re = tbl["Potential"]["Morse"]["re"].value_or(re);
      cutoff = tbl["Potential"]["Morse"]["cutoff"].value_or(cutoff);
    }
  };

private:
  void morse(double r, double &energy, double &force);
  double _De;
  double _a;
  double _re;
  double _cutoff;
  double _energyCutoff;

public:
  Morse(const Morse::Params &p_a)
      : _De{p_a.De},
        _a{p_a.a},
        _re{p_a.re},
        _cutoff{p_a.cutoff} {}
  void forceImpl(const ForceInput &, ForceOut *) override final;
  void setParameters(const Morse::Params &p_a);
};

} // namespace eonc
