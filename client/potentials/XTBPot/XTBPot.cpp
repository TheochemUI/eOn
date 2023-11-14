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
using forcefields::unit_system::HARTREE;
using forcefields::unit_system::BOHR;

// pointer to number of atoms, pointer to array of positions
// pointer to array of forces, pointer to internal energy
// address to supercell size
void XTBPot::force(long N, const double *R, const int *atomicNrs, double *F,
                   double *U, double *variance, const double *box) {
  variance = nullptr;
  int intN = static_cast<int>(N);
  const bool periodicity[3]{false, false, false};

  // Allocate memory for converted positions
  double R_bohr[3 * N];

  // Convert positions from Angstroms to Bohrs
  for (long i = 0; i < 3 * N; ++i) {
    R_bohr[i] = R[i] / BOHR;
  }

  xtb_TMolecule mol = xtb_newMolecule(env, &intN, atomicNrs, R_bohr, nullptr,
                                      nullptr, nullptr, nullptr);
  if (!mol) {
    throw std::runtime_error("Failed to create xtb molecule");
  }

  // Load a specific GFN-xTB calculator
  // Ordered from lowest accuracy to highest
  // xtb_loadGFNFF(env, mol, calc, NULL);
  // xtb_loadGFN0xTB(env, mol, calc, NULL);
  // xtb_loadGFN1xTB(env, mol, calc, NULL);
  xtb_loadGFN2xTB(env, mol, calc, NULL);
  xtb_setAccuracy(env, calc, 1);
  xtb_setElectronicTemp(env, calc, 300.0);
  xtb_setMaxIter(env, calc, 250);

  // Calculate
  xtb_TResults res = xtb_newResults();
  if (!res) {
    xtb_delMolecule(&mol);
    throw std::runtime_error("Failed to create xtb results");
  }

  xtb_singlepoint(env, mol, calc, res);

  // Extract energy
  xtb_getEnergy(env, res, U);

  // Extract forces
  xtb_getGradient(env, res, F);

  // Convert back to angstrom and eV based units
  for (int i = 0; i < N; i++) {
    F[3 * i] *= -1*(HARTREE / BOHR);
    F[3 * i + 1] *= -1*(HARTREE / BOHR);
    F[3 * i + 2] *= -1*(HARTREE / BOHR);
  }
  *U *= HARTREE;

  // Clean up molecule and results objects
  xtb_delResults(&res);
  xtb_delMolecule(&mol);
}
