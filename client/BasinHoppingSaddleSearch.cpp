#include "BasinHoppingSaddleSearch.h"
#include "Dimer.h"
#include "ImprovedDimer.h"
#include "Lanczos.h"
#include "Log.h"
#include "LowestEigenmode.h"
#include "MinModeSaddleSearch.h"
#include "NudgedElasticBand.h"
#include <stdio.h>

int BasinHoppingSaddleSearch::run(void) {
  // minimize "saddle"
  saddle->relax(false, true, false, "displacementmin");
  product = std::make_shared<Matter>(pot, params);
  *product = *saddle;
  // accept or reject based on boltzman exp(-de/(kB*params->temperature))
  double eproduct, ereactant, de;
  eproduct = product->getPotentialEnergy();
  ereactant = reactant->getPotentialEnergy();
  de = eproduct - ereactant;
  double kB = params->kB;
  double Temperature = params->temperature;
  double arg = -de / (kB * Temperature);
  double p = exp(arg);
  double r = helper_functions::random();
  if (ereactant < eproduct) {
    if (r > p) { // reject
      return 1;
    }
  }
  // NEB reactant to minimized "saddle"
  NudgedElasticBand neb(reactant, product, params, pot);
  neb.image[0]->matter2con("neb_initial_band.con", false);
  for (int j = 1; j < neb.images; j++) {
    neb.image[j]->matter2con("neb_initial_band", true);
  }
  neb.compute();
  // pick the maximum energy image along the band
  double Emax = -1e100;
  int HighestImage = 0;

  for (int i = 1; i < neb.images; i++) {
    double Etest = neb.image[i]->getPotentialEnergy();
    printf("i: %i Etest: %f \n", i, Etest);
    if (Etest > Emax) {
      Emax = Etest;
      HighestImage = i;
    }
  }
  // do dimer
  // Calculate initial direction
  AtomMatrix r_1 = neb.image[HighestImage - 1]->getPositions();
  AtomMatrix r_2 = neb.image[HighestImage]->getPositions();
  AtomMatrix r_3 = neb.image[HighestImage + 1]->getPositions();
  AtomMatrix direction = (r_3 - r_1) / 2;
  MinModeSaddleSearch dim(neb.image[HighestImage], direction.normalized(),
                          ereactant, params, pot);
  dim.run();
  *saddle = *neb.image[HighestImage];
  eigenvalue = dim.getEigenvalue();
  eigenvector = dim.getEigenvector();
  return 0;
}

double BasinHoppingSaddleSearch::getEigenvalue() { return eigenvalue; }

AtomMatrix BasinHoppingSaddleSearch::getEigenvector() { return eigenvector; }
