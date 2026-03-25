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

#include <cstring>
#include <vector>

void EffectiveMediumTheory::cleanMemory() {
  delete EMTObj;
  EMTObj = nullptr;
  delete SuperCellObj;
  SuperCellObj = nullptr;
  delete AtomsObj;
  AtomsObj = nullptr;
  delete EMTParameterObj;
  EMTParameterObj = nullptr;
}

// pointer to number of atoms, pointer to array of positions
// pointer to array of forces, pointer to internal energy
// adress to supercell size
void EffectiveMediumTheory::force(long N, const double *R, const int *atomicNrs,
                                  double *F, double *U, double *variance,
                                  const double *box) {
  variance = nullptr;

  // Copy positions (Asap may modify them internally)
  std::vector<double> pos(R, R + 3 * N);

  // Reinitialize if atom count changed
  if (numberOfAtoms != N) {
    numberOfAtoms = N;
    cleanMemory();

    Vec tempBasis[3] = {
        Vec(box[0], box[1], box[2]),
        Vec(box[3], box[4], box[5]),
        Vec(box[6], box[7], box[8]),
    };

    SuperCellObj = new SuperCell(tempBasis, periodicity);
    AtomsObj = new Atoms(reinterpret_cast<Vec *>(pos.data()), N, SuperCellObj);

    std::vector<int> atomicNrsTemp(atomicNrs, atomicNrs + N);
    AtomsObj->SetAtomicNumbers(atomicNrsTemp.data());

    if (emtRasmussen) {
      EMTParameterObj = new EMTRasmussenParameterProvider();
      EMTObj = new EMT(EMTParameterObj);
    } else {
      EMTObj = new EMT(nullptr);
    }
    AtomsObj->SetCalculator(EMTObj);
  }

  AtomsObj->SetCartesianPositions(reinterpret_cast<Vec *>(pos.data()));

  // Update the box
  Vec tempBasis[3] = {
      Vec(box[0], box[1], box[2]),
      Vec(box[3], box[4], box[5]),
      Vec(box[6], box[7], box[8]),
  };
  AtomsObj->SetUnitCell(tempBasis, true);

  *U = EMTObj->GetPotentialEnergy();

  const Vec *tempF = EMTObj->GetCartesianForces();
  std::memcpy(F, tempF, N * sizeof(Vec));
}
