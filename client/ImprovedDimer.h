#ifndef IMPROVEDDIMER_H
#define IMPROVEDDIMER_H

#include "Eigen.h"
#include "LowestEigenmode.h"
#include "Matter.h"
#include "Parameters.h"
#include <vector>

// dimer method to find the lowest curvature mode
class ImprovedDimer : public LowestEigenmode {

private:
  shared_ptr<spdlog::logger> log;

public:
  // Optimization for the dimer
  //    static const string OPT_SD;
  //    static const string OPT_CG;
  //    static const string OPT_LBFGS;
  static const char OPT_SD[];
  static const char OPT_CG[];
  static const char OPT_LBFGS[];

  ImprovedDimer(std::shared_ptr<Matter> matter,
                std::shared_ptr<Parameters> params,
                std::shared_ptr<Potential> pot);
  ~ImprovedDimer() = default;

  void compute(std::shared_ptr<Matter> matter, AtomMatrix initialDirection);
  double getEigenvalue();
  AtomMatrix getEigenvector();

  std::shared_ptr<Matter> x0; // Center image
  std::shared_ptr<Matter> x1; // Forward image
  VectorXd tau;               // Dimer direction
  VectorXd theta;             // Dimer rotation direction
  VectorXd F_R;               // Dimer rotational force
  double C_tau;               // Curvature along tau

  // parameters used for conjugate gradients
  VectorXd F_R_Old;
  VectorXd thetaOld;
  double a, b, gamma;
  bool init_cg;

  // variables for LBFGS
  std::vector<VectorXd> s, y;
  std::vector<double> rho;
  bool init_lbfgs;
  VectorXd rPrev;

  std::vector<VectorXd> gradients;
  std::vector<VectorXd> positions;
};

#endif
