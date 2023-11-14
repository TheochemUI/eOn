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

// Conversion factor from Angstrom to Bohr
const double angstromToBohr = 1.88973;
const double hartreeToEV = 27.2114;

// pointer to number of atoms, pointer to array of positions
// pointer to array of forces, pointer to internal energy
// address to supercell size
void XTBPot::force(long N, const double *R, const int *atomicNrs, double *F,
                   double *U, double *variance, const double *box) {
  variance = nullptr;
  *U = 0;
  for (int i = 0; i < N; i++) {
    F[3 * i] = 0;
    F[3 * i + 1] = 0;
    F[3 * i + 2] = 0;
  }

  int intN = static_cast<int>(N);
  const bool periodic[3]  {false, false, false};

  // Allocate memory for converted positions
  double* R_bohr = new double[3 * N];

  // Convert positions from Angstroms to Bohrs
  for (long i = 0; i < 3 * N; ++i) {
    R_bohr[i] = R[i] * angstromToBohr;
  }

  xtb_TMolecule mol = xtb_newMolecule(env, &intN, atomicNrs, R_bohr, nullptr,
                                      nullptr, box, periodic);
  if (!mol) {
    delete[] R_bohr;
    throw std::runtime_error("Failed to create xtb molecule");
  }

  // Load a specific GFN-xTB calculator
  xtb_loadGFN1xTB(env, mol, calc, nullptr);

  xtb_TResults res = xtb_newResults();
  if (!res) {
    xtb_delMolecule(&mol);
    delete[] R_bohr;
    throw std::runtime_error("Failed to create xtb results");
  }

  xtb_singlepoint(env, mol, calc, res);

  // Extract energy
  xtb_getEnergy(env, res, U);

  // Extract forces
  // std::vector<double> forces(N * 3);
  xtb_getGradient(env, res, F);

  // Convert back to angstrom and eV
  for (int i = 0; i < N; i++) {
    F[3 * i] /= angstromToBohr;
    F[3 * i + 1] /= angstromToBohr;
    F[3 * i + 2] /= angstromToBohr;
  }
  *U *= hartreeToEV;


  // Copy forces to F
  // std::copy(forces.begin(), forces.end(), F);

  // Clean up molecule and results objects
  xtb_delResults(&res);
  xtb_delMolecule(&mol);
  delete[] R_bohr;
}
