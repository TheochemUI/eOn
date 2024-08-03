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

#include "ExtPot.h"
#include <iomanip>
#include <stdio.h>
#include <unistd.h>
namespace eonc {

void ExtPot::forceImpl(const ForceInput &fip, ForceOut *efvd) {
#ifdef EON_CHECKS
  eonc::pot::checkParams(fip);
  eonc::pot::zeroForceOut(fip.nAtoms, efvd);
#endif
  const long int N = fip.nAtoms;
  passToSystem(N, fip.pos, fip.atmnrs, fip.box);
  int result = std::system(eon_extpot_path.c_str());
  if (result != 0) {
    throw std::runtime_error("Error executing external potential system call.");
  }
  receiveFromSystem(N, efvd->F, &efvd->energy);
  return;
}

void ExtPot::passToSystem(long N, const double *R, const size_t *atomicNrs,
                          const double *box) {
  std::ofstream out("from_eon_to_extpot");
  if (!out) {
    throw std::runtime_error(
        "Could not open 'from_eon_to_extpot' for writing.");
  }

  for (int i = 0; i < 3; ++i) {
    out << std::fixed << std::setprecision(19) << box[i * 3 + 0] << "\t"
        << box[i * 3 + 1] << "\t" << box[i * 3 + 2] << "\n";
  }

  for (long i = 0; i < N; ++i) {
    out << atomicNrs[i] << "\t" << R[i * 3 + 0] << "\t" << R[i * 3 + 1] << "\t"
        << R[i * 3 + 2] << "\n";
  }
}

void ExtPot::receiveFromSystem(long N, double *F, double *U) {
  std::ifstream in("from_extpot_to_eon");
  if (!in) {
    throw std::runtime_error(
        "Could not open 'from_extpot_to_eon' for reading.");
  }

  in >> *U;

  for (long i = 0; i < N; ++i) {
    in >> F[i * 3 + 0] >> F[i * 3 + 1] >> F[i * 3 + 2];
  }
}

} // namespace eonc
