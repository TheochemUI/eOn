//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include "GPRPotential.h"
#include "../../subprojects/gprdimer/gpr/auxiliary/AdditionalFunctionality.h"
#include "../../subprojects/gprdimer/structures/Structures.h"

namespace {

const char *elementArray[] = {
    "Unknown", "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na",
    "Mg",      "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",
    "Cr",      "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br",
    "Kr",      "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag",
    "Cd",      "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr",
    "Nd",      "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
    "Hf",      "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi",
    "Po",      "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",  NULL};

// guess the atom type from the atomic mass,
std::string mass2atom(double atomicmass) {
  return elementArray[int(atomicmass + .5)];
}

int symbol2atomicNumber(char const *symbol) {
  int i = 0;

  while (elementArray[i] != NULL) {
    if (strcmp(symbol, elementArray[i]) == 0) {
      return i;
    }
    i++;
  }
  // invalid symbol
  return -1;
}

char const *atomicNumber2symbol(int n) { return elementArray[n]; }
} // namespace

GPRPotential::GPRPotential(Parameters *p) { gpr_model = nullptr; }

void GPRPotential::registerGPRObject(
    gpr::GaussianProcessRegression *_gpr_model) {
  gpr_model = _gpr_model;
}

void GPRPotential::initialize(void) {}

void GPRPotential::cleanMemory(void) {}

// pointer to number of atoms, pointer to array of positions
// pointer to array of forces, pointer to internal energy
// adress to supercell size
void GPRPotential::force(long N, const double *R, const int *atomicNrs,
                         double *F, double *U, double *variance,
                         const double *box) {
  variance = nullptr;
  gpr::Observation observation;

  // Copy R points. Note, R should correspond to the moving atoms only.
  observation.R.resize(1, N * 3);
  for (int i = 0; i < N; i++) {
    observation.R.set(i, {R[3 * i], R[3 * i + 1], R[3 * i + 2]});
  }

  // Note, the following functions should be called before calling for
  // gpr_model->calculatePotential() gpr_model->decomposeCovarianceMatrix(R,
  // ind) - takes covariance matrix and vector of repetitive indices
  // gpr_model->calculateMeanPrediction() - takes a vector of combined energy
  // and force gpr_model->calculatePosteriorMeanPrediction() - no arguments
  gpr_model->calculatePotential(observation);

  for (int i = 0; i < N; i++) {
    F[3 * i] = observation.G[3 * i];
    F[3 * i + 1] = observation.G[3 * i + 1];
    F[3 * i + 2] = observation.G[3 * i + 2];
  }

  // FIXME: Test conversion, E should only have one element here
  *U = observation.E[0];
}
