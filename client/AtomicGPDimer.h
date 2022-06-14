#ifndef ATOMICGPDIMER_H
#define ATOMICGPDIMER_H

#include "Eigen.h"
#include "LowestEigenmode.h"
#include "Matter.h"
#include "Parameters.h"
#include <vector>

// Defined in the gprdimer target
#include "subprojects/gprdimer/gpr/saddle_point/AtomicDimer.h"
#include "subprojects/gprdimer/gpr/Enums.h"
#include "subprojects/gprdimer/structures/Structures.h"
#include "subprojects/gprdimer/gpr/auxiliary/ProblemSetUp.h"
#include "subprojects/gprdimer/gpr/covariance_functions/ConstantCF.h"
#include "subprojects/gprdimer/gpr/covariance_functions/SexpatCF.h"
#include "subprojects/gprdimer/gpr/ml/GaussianProcessRegression.h"
#include "subprojects/gprdimer/managers/io/FileManager.h"
#include "subprojects/gprdimer/structures/Structures.h"

// dimer method to find the lowest curvature mode
class AtomicGPDimer : public LowestEigenmode {

public:
  // Optimization for the dimer
  static const char OPT_SCG[];
  static const char OPT_LBFGS[];

  AtomicGPDimer(Matter *matter, Parameters *parameters);
  ~AtomicGPDimer();

  void compute(Matter *matter, AtomMatrix initialDirectionMatrix);
  double getEigenvalue();
  AtomMatrix getEigenvector();

private:
  Matter *matterCenter;       // initial center of the dimer
  AtomMatrix direction;       // direction along the dimer
  AtomMatrix rotationalPlane; // direction normal to the plane of dimer rotation
  Parameters *parameters;

  gpr::InputParameters p;
  atmd::AtomicDimer atomic_dimer;
  aux::ProblemSetUp problem_setup;
  gpr::AtomsConfiguration atoms_config;
  gpr::Observation init_observations, init_middle_point;
  gpr::Coord orient_init, R_init;
};

#endif
