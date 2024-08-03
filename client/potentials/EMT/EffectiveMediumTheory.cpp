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

#include "EffectiveMediumTheory.h"
#include <string.h>
namespace eonc {
void EffectiveMediumTheory::cleanMemory() {
  EMTObj.reset();
  SuperCellObj.reset();
  AtomsObj.reset();
  EMTParameterObj.reset();
}

void EffectiveMediumTheory::initialize(const ForceInput &fip) {
  const size_t N = fip.nAtoms;
  numberOfAtoms = N;

  pos.resize(N);
  atomicNrsTemp.resize(N);
  for (size_t i = 0; i < N; ++i) {
    pos[i] = Vec(fip.pos[3 * i], fip.pos[3 * i + 1], fip.pos[3 * i + 2]);
    atomicNrsTemp[i] = static_cast<int>(fip.atmnrs[i]);
  }

  tempBasis = {Vec(fip.box[0], fip.box[1], fip.box[2]),
               Vec(fip.box[3], fip.box[4], fip.box[5]),
               Vec(fip.box[6], fip.box[7], fip.box[8])};

  SuperCellObj = std::make_unique<SuperCell>(tempBasis.data(), periodicity);
  AtomsObj = std::make_unique<Atoms>(pos.data(), static_cast<int>(N),
                                     SuperCellObj.get());
  AtomsObj->SetAtomicNumbers(atomicNrsTemp.data());

  if (useEMTRasmussen) {
    EMTParameterObj = std::make_unique<EMTRasmussenParameterProvider>();
    EMTObj = std::make_unique<EMT>(EMTParameterObj.get());
  } else {
    EMTObj = std::make_unique<EMT>(nullptr);
  }

  AtomsObj->SetCalculator(EMTObj.get());
}

void EffectiveMediumTheory::forceImpl(const ForceInput &fip, ForceOut *efvd) {
#ifdef EON_CHECKS
  eonc::pot::checkParams(fip);
  eonc::pot::zeroForceOut(fip.nAtoms, efvd);
#endif
  const size_t N = fip.nAtoms;

  if (numberOfAtoms != N) {
    cleanMemory();
    initialize(fip);
  } else {
    pos.resize(N);
    for (size_t i = 0; i < N; ++i) {
      pos[i] = Vec(fip.pos[3 * i], fip.pos[3 * i + 1], fip.pos[3 * i + 2]);
    }
    AtomsObj->SetCartesianPositions(pos.data());
  }

  tempBasis = {Vec(fip.box[0], fip.box[1], fip.box[2]),
               Vec(fip.box[3], fip.box[4], fip.box[5]),
               Vec(fip.box[6], fip.box[7], fip.box[8])};
  AtomsObj->SetUnitCell(tempBasis.data(), true);

  efvd->energy = EMTObj->GetPotentialEnergy();

  const auto &tempF = EMTObj->GetCartesianForces();
  std::memcpy(efvd->F, tempF, N * sizeof(Vec));
}

} // namespace eonc
