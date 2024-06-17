#ifndef ATOMICGPDIMER_H
#define ATOMICGPDIMER_H

#include "Eigen.h"
#include "LowestEigenmode.h"
#include "Matter.h"
#include "Parameters.h"
#include <vector>

// Defined in the gprdimer target
#include "subprojects/gprdimer/gpr/AtomicDimer.h"
#include "subprojects/gprdimer/gpr/Enums.h"
#include "subprojects/gprdimer/gpr/auxiliary/ProblemSetUp.h"
#include "subprojects/gprdimer/gpr/covariance_functions/ConstantCF.h"
#include "subprojects/gprdimer/gpr/covariance_functions/SexpatCF.h"
#include "subprojects/gprdimer/gpr/ml/GaussianProcessRegression.h"
#include "subprojects/gprdimer/managers/io/FileManager.h"

// dimer method to find the lowest curvature mode
class AtomicGPDimer : public LowestEigenmode {

public:
  // Optimization for the dimer
  static const char OPT_SCG[];
  static const char OPT_LBFGS[];

  AtomicGPDimer(std::shared_ptr<Matter> matter,
                std::shared_ptr<Parameters> params,
                std::shared_ptr<Potential> pot);
  ~AtomicGPDimer() = default;

  void compute(std::shared_ptr<Matter> matter, AtomMatrix initialDirection);
  double getEigenvalue();
  AtomMatrix getEigenvector();

private:
  std::shared_ptr<Matter> matterCenter; // initial center of the dimer
  AtomMatrix direction;                 // direction along the dimer
  AtomMatrix rotationalPlane; // direction normal to the plane of dimer rotation

  gpr::InputParameters p;
  atmd::AtomicDimer atomic_dimer;
  aux::ProblemSetUp problem_setup;
  gpr::AtomsConfiguration atoms_config;
  gpr::Observation init_observations, init_middle_point;
  gpr::Coord orient_init, R_init;
};

#endif
