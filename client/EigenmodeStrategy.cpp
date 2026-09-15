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

#include "eon/EigenmodeStrategy.h"
#include "eon/Davidson.h"
#include "eon/Dimer.h"
#include "eon/ImprovedDimer.h"
#include "eon/Lanczos.h"
#include "eon/Parameters.h"
#include <stdexcept>

#ifdef WITH_GPRD
#include "eon/AtomicGPDimer.h"
#endif

namespace eonc {

std::shared_ptr<LowestEigenmode>
buildEigenmodeStrategy(std::shared_ptr<Matter> matter, const Parameters &params,
                       std::shared_ptr<Potential> pot) {
  if (params.saddle_search_options.minmode_method ==
      LowestEigenmode::MINMODE_DIMER) {
    if (params.dimer_options.improved) {
      return std::make_shared<ImprovedDimer>(matter, params, pot);
    }
    return std::make_shared<Dimer>(matter, params, pot);
  }
  if (params.saddle_search_options.minmode_method ==
      LowestEigenmode::MINMODE_LANCZOS) {
    return std::make_shared<Lanczos>(matter, params, pot);
  }
  if (params.saddle_search_options.minmode_method ==
      LowestEigenmode::MINMODE_DAVIDSON) {
    return std::make_shared<Davidson>(matter, params, pot);
  }
  if (params.saddle_search_options.minmode_method ==
      LowestEigenmode::MINMODE_GPRDIMER) {
#ifdef WITH_GPRD
    return std::make_shared<AtomicGPDimer>(matter, params, pot);
#else
    throw std::runtime_error(
        "min_mode_method=gprdimer requires -Dwith_gprd=true (WITH_GPRD)");
#endif
  }
  return std::make_shared<ImprovedDimer>(matter, params, pot);
}

} // namespace eonc
