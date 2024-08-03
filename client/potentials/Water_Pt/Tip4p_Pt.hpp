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
Wrapper for Eon
@author Jean-Claude C. Berthet
@date 2007
University of Iceland
*/

#pragma once
#include "../../Potential.h"
#include "zhu_philpott.hpp"
namespace eonc {
class Tip4p_Pt : public Potential<Tip4p_Pt>,
                 private forcefields::ZhuPhilpott<> {
public:
  Tip4p_Pt()
      : // TODO(rg): Expose these
        forcefields::ZhuPhilpott<>(8.5, 1.0) {}
  void forceImpl(const ForceInput &, ForceOut *) override final;
};
} // namespace eonc
