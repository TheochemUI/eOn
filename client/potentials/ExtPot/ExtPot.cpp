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

#include "eon/potentials/ExtPot/ExtPot.h"
#include "eon/potentials/ExternalCommand.h"

#include <format>
#include <fstream>
#include <stdexcept>
#include <string>

namespace {

// The exchange files live in the current working directory and carry fixed
// names, so two clients running in one directory read each other's numbers.
constexpr const char *kToExtPot = "from_eon_to_extpot";
constexpr const char *kFromExtPot = "from_extpot_to_eon";

} // namespace

void ExtPot::cleanMemory(void) { return; }

ExtPot::~ExtPot() { cleanMemory(); }

void ExtPot::force(long N, const double *R, const int *atomicNrs, double *F,
                   double *U, double *variance, const double *box) {
  variance = nullptr;
  if (eon_extpot_path.empty()) {
    throw std::runtime_error(
        "ExtPot needs potential_options.extPotPath to name a program to run");
  }
  passToSystem(N, R, atomicNrs, box);
  eonc::pot::runOrThrow(eon_extpot_path);
  recieveFromSystem(N, F, U);
  return;
}

void ExtPot::passToSystem(long N, const double *R, const int *atomicNrs,
                          const double *box)
// 'positions' of all particles and box
{
  std::ofstream out(kToExtPot, std::ios::trunc);
  if (!out) {
    throw std::runtime_error(
        std::format("Could not open {} for writing", kToExtPot));
  }

  for (int i = 0; i < 3; i++) {
    out << std::format("{:.19f}\t{:.19f}\t{:.19f}\n", box[i * 3 + 0],
                       box[i * 3 + 1], box[i * 3 + 2]);
  }

  for (long i = 0; i < N; i++) {
    out << std::format("{}\t{:.19f}\t{:.19f}\t{:.19f}\n", atomicNrs[i],
                       R[i * 3 + 0], R[i * 3 + 1], R[i * 3 + 2]);
  }

  out.close();
  if (!out) {
    throw std::runtime_error(
        std::format("Could not write the structure to {}", kToExtPot));
  }
  return;
}

void ExtPot::recieveFromSystem(long N, double *F, double *U)
// first line must be the total 'energy', the following lines should be the
// 'forces'
{
  std::ifstream in(kFromExtPot);
  if (!in) {
    throw std::runtime_error(std::format(
        "Could not open {}; the external program left no result", kFromExtPot));
  }

  if (!(in >> *U)) {
    throw std::runtime_error(
        std::format("Could not read the energy from {}", kFromExtPot));
  }

  for (long i = 0; i < N; i++) {
    if (!(in >> F[i * 3 + 0] >> F[i * 3 + 1] >> F[i * 3 + 2])) {
      throw std::runtime_error(
          std::format("{} holds forces for {} atoms, expected {}", kFromExtPot,
                      i, N));
    }
  }
  return;
}
