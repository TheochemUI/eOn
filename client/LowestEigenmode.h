#ifndef LOWESTEIGENMODE_H
#define LOWESTEIGENMODE_H

#include "Eigen.h"
#include "Matter.h"
#include "Parameters.h"

/* Define the interface for the lowest eigenvalue determination algorithm */
class LowestEigenmode {
protected:
  enum class LEMStatus {
    GOOD = 0,
    INIT = 1,
    BAD_MAX_ITERATIONS = 2,
    RUNNING,
    MAX_UNCERTAINITY
  };
  // make const
  std::shared_ptr<Potential> pot;
  std::shared_ptr<Parameters> params;
  LowestEigenmode::LEMStatus status;

public:
  // stats information
  long totalForceCalls;
  double statsTorque;
  double statsCurvature;
  double statsAngle;
  long statsRotations;
  long totalIterations; // Only set by the gpr dimer
  static const char MINMODE_DIMER[];
  static const char MINMODE_GPRDIMER[];
  static const char MINMODE_LANCZOS[];

  LowestEigenmode(std::shared_ptr<Potential> potPassed,
                  std::shared_ptr<Parameters> parameters)
      : pot{potPassed},
        params{parameters} {}
  virtual ~LowestEigenmode() = default;

  // void virtual initialize(Matter const *matter, AtomMatrix displacement) = 0;
  virtual LowestEigenmode::LEMStatus compute(std::shared_ptr<Matter> matter,
                                             AtomMatrix direction) = 0;
  LowestEigenmode::LEMStatus getStatus() { return this->status; };
  virtual double getEigenvalue() = 0;
  virtual AtomMatrix getEigenvector() = 0;
};

#endif
