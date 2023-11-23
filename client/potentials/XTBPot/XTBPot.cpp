//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

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

  // Make molecule
  xtb_TMolecule mol = xtb_newMolecule(env, &intN, atomicNrs, R_bohr, nullptr,
                                      nullptr, box_bohr, periodicity);

  // Setup parameters
  if (xtb_paramset == GFNMethod::GFNFF) {
    xtb_loadGFNFF(env, mol, calc, nullptr);
  } else if (xtb_paramset == GFNMethod::GFN0xTB) {
    xtb_loadGFN0xTB(env, mol, calc, nullptr);
  } else if (xtb_paramset == GFNMethod::GFN1xTB) {
    xtb_loadGFN1xTB(env, mol, calc, nullptr);
  } else if (xtb_paramset == GFNMethod::GFN2xTB) {
    xtb_loadGFN2xTB(env, mol, calc, nullptr);
  } else {
    throw std::runtime_error("Parameter set for XTB must be one of GFNFF, "
                             "GFN0xTB, GFN1xTB or GFN2xTB.\n");
  }
  xtb_setAccuracy(env, calc, xtb_acc);
  xtb_setElectronicTemp(env, calc, xtb_electronic_temperature);
  xtb_setMaxIter(env, calc, xtb_max_iter);

  // Calculate
  xtb_TResults res = xtb_newResults();
  xtb_singlepoint(env, mol, calc, res);
  counter++;

  // Extract results
  xtb_getEnergy(env, res, U);
  xtb_getGradient(env, res, F);

  // Convert back to angstrom and eV based units
  for (int idx = 0; idx < N; idx++) {
    F[3 * idx] *= -1 * (HARTREE / BOHR);
    F[3 * idx + 1] *= -1 * (HARTREE / BOHR);
    F[3 * idx + 2] *= -1 * (HARTREE / BOHR);
  }
  *U *= HARTREE;

  // Clean up results
  xtb_delResults(&res);
  xtb_delMolecule(&mol);
}
