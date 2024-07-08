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

// General Functions
void EffectiveMediumTheory::cleanMemory(void) {
  if (EMTObj != 0) {
    delete EMTObj;
    EMTObj = 0;
  }
  if (SuperCellObj != 0) {
    delete SuperCellObj;
    SuperCellObj = 0;
  }
  if (AtomsObj != 0) {
    delete AtomsObj;
    AtomsObj = 0;
  }
  if (EMTParameterObj != 0) {
    delete EMTParameterObj;
    EMTParameterObj = 0;
  }
  return;
}

// pointer to number of atoms, pointer to array of positions
// pointer to array of forces, pointer to internal energy
// adress to supercell size
void EffectiveMediumTheory::force(long N, const double *R, const int *atomicNrs,
                                  double *F, double *U, double *variance,
                                  const double *box) {
  variance = nullptr;
  int i, j;
  double *pos;

  pos = new double[3 * N];

  for (i = 0; i < 3 * N; i++)
    pos[i] = R[i];

  // an atom has been deposited
  if (numberOfAtoms != N) {
    numberOfAtoms = N;
    cleanMemory();
    int *atomicNrsTemp;
    atomicNrsTemp = new int[N];

    // create new atoms / emt potential with N atoms
    Vec tempBasisX(box[0], box[1], box[2]);
    Vec tempBasisY(box[3], box[4], box[5]);
    Vec tempBasisZ(box[6], box[7], box[8]);
    Vec tempBasis[3];
    tempBasis[0] = tempBasisX;
    tempBasis[1] = tempBasisY;
    tempBasis[2] = tempBasisZ;

    SuperCellObj = new SuperCell(tempBasis, periodicity);
    AtomsObj = new Atoms((Vec *)pos, N, SuperCellObj);

    for (j = 0; j < N; j++)
      atomicNrsTemp[j] = int(atomicNrs[j]);

    AtomsObj->SetAtomicNumbers(atomicNrsTemp);

    if (m_params->EMTRasmussen) {
      EMTParameterObj = new EMTRasmussenParameterProvider();
      EMTObj = new EMT(EMTParameterObj);
    } else
      EMTObj = new EMT(NULL);

    AtomsObj->SetCalculator(EMTObj);
    delete[] atomicNrsTemp;
  }
  AtomsObj->SetCartesianPositions((Vec *)pos);
  // update the box
  Vec tempBasisX(box[0], box[1], box[2]);
  Vec tempBasisY(box[3], box[4], box[5]);
  Vec tempBasisZ(box[6], box[7], box[8]);
  Vec tempBasis[3];
  tempBasis[0] = tempBasisX;
  tempBasis[1] = tempBasisY;
  tempBasis[2] = tempBasisZ;
  AtomsObj->SetUnitCell(tempBasis, true);

  *U = EMTObj->GetPotentialEnergy();

  // converts data from EMT to suite EON
  const Vec *tempF = EMTObj->GetCartesianForces();
  memcpy(F, tempF, N * sizeof(Vec));

  delete[] pos;
  return;
}
