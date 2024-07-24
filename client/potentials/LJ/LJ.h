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
  friend class LJCluster;

public:
  struct Params final {
    double u0{1.0};
    double cutoff_R{15.0};
    double psi{1.0};
  };

private:
  double _u0;
  double _cutoff_R;
  double _psi;
  double _cutoff_U{0};

  double calc_cutoffU(const LJ::Params &);

public:
  LJ(const LJ::Params &p_a)
      : _u0{p_a.u0},
        _cutoff_R{p_a.cutoff_R},
        _psi{p_a.psi},
        _cutoff_U{calc_cutoffU(p_a)} {}

  void forceImpl(const ForceInput &, ForceOut *) override final;

  void setParameters(const LJ::Params &);
};

} // namespace eonc
