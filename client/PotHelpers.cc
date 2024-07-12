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
#include "PotHelpers.hpp"
#include <stdexcept>

// These are all only run if EON_CHECKS is set, i.e. not in release mode
namespace eonc::pot {

void zeroForceOut(const size_t &nAtoms, ForceOut *efvd) {
  efvd->energy = 0;
  efvd->variance = 0;
  for (size_t idx{0}; idx < nAtoms; idx++) {
    efvd->F[3 * idx] = 0;
    efvd->F[3 * idx + 1] = 0;
    efvd->F[3 * idx + 2] = 0;
  }
};

void checkParams(const ForceInput &params) {
  // TODO(rg):: Can't the box be in negative cartesian sapce?
  if (!(params.box[0] > 0) and (params.box[4] > 0) and (params.box[8] > 0)) {
    // For Morse, might be specific to it only
    throw std::runtime_error("Can't work with non-zero box in force call");
  }
  if (params.nAtoms <= 0) {
    throw std::runtime_error("Can't work with zero atoms in force call");
  }
}

} // namespace eonc::pot
