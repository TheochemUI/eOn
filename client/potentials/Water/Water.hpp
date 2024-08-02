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
/** @file
Wrapper for eOn
@author Jean-Claude C. Berthet
@date 2007
University of Iceland
*/

#pragma once
#include "../../Potential.h"
#include "spce_ccl.hpp"
#include "tip4p_ccl.hpp"
namespace eonc {
class Tip4p : public Potential<Tip4p>, private forcefields::Tip4p {
public:
  Tip4p()
      : // TODO(rg): Expose these like LJ
        forcefields::Tip4p(8.5, 1.0) {};
  void forceImpl(const ForceInput &, ForceOut *) override final;
};

class SpceCcl : public Potential<SpceCcl>, private forcefields::SpceCcl {
public:
  SpceCcl()
      : forcefields::SpceCcl(8.5, 1.0) {}
  void forceImpl(const ForceInput &, ForceOut *) override final;
};

} // namespace eonc
