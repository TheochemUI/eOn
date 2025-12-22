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

#include "XTBPot.h"
#include <cstddef>

// Conversion factors
// const double angstromToBohr = 1.8897261349925714;
// const double hartreeToEV = 27.21138386;
// const double hartreeBohr_to_eVA = 14.399645472115932;
using forcefields::unit_system::BOHR;
using forcefields::unit_system::HARTREE;

// pointer to number of atoms, pointer to array of positions
// pointer to array of forces, pointer to internal energy
// address to supercell size
void XTBPot::force(long N, const double *R, const int *atomicNrs, double *F,
                   double *U, double *variance, const double *box) {
  variance = nullptr;
  int intN = static_cast<int>(N);
  // TODO: Periodicity shouldn't crash
  const bool periodicity[3]{false, false, false};
  double box_bohr[3 * 3];

  // Allocate memory for converted positions
  double R_bohr[3 * N];

  // Convert positions from Angstrom to Bohr
  for (long idx = 0; idx < 3 * N; ++idx) {
    R_bohr[idx] = R[idx] / BOHR;
  }
  for (long idx = 0; idx < 9; ++idx) {
    box_bohr[idx] = box[idx] / BOHR;
  }

  if (!initialized) {
    // First call: Create the molecule and load the Hamiltonian
    mol = xtb_newMolecule(env, &intN, atomicNrs, R_bohr, &total_charge, &uhf,
                          box_bohr, periodicity);

    switch (xtb_paramset) {
    case GFNMethod::GFNFF:
      xtb_loadGFNFF(env, mol, calc, nullptr);
      break;
    case GFNMethod::GFN0xTB:
      xtb_loadGFN0xTB(env, mol, calc, nullptr);
      break;
    case GFNMethod::GFN1xTB:
      xtb_loadGFN1xTB(env, mol, calc, nullptr);
      break;
    case GFNMethod::GFN2xTB:
      xtb_loadGFN2xTB(env, mol, calc, nullptr);
      break;
    }

    xtb_setAccuracy(env, calc, xtb_acc);
    xtb_setElectronicTemp(env, calc, xtb_electronic_temperature);
    xtb_setMaxIter(env, calc, xtb_max_iter);
    initialized = true;
  } else {
    // Subsequent calls: Only update coordinates and lattice
    xtb_updateMolecule(env, mol, R_bohr, box_bohr);
  }

  xtb_singlepoint(env, mol, calc, res);

  // Check for SCF convergence or internal xTB errors
  if (xtb_checkEnvironment(env) != 0) {
    char err_msg[512];
    xtb_getError(env, err_msg, nullptr);
    throw std::runtime_error(std::string("xTB Error: ") + err_msg);
  }

  xtb_getEnergy(env, res, U);
  xtb_getGradient(env, res, F);

  // Convert Hartree/Bohr to eV/Angstrom
  for (long i = 0; i < 3 * N; ++i) {
    F[i] *= -1.0 * (HARTREE / BOHR);
  }
  *U *= HARTREE;
  counter++;
}
