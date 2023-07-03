#ifndef OBJECTIVE_FUNCTION_H
#define OBJECTIVE_FUNCTION_H

#include "Dimer.h"
#include "ImprovedDimer.h"
#include "Lanczos.h"
#include "LowestEigenmode.h"
#include "Matter.h"

#ifdef WITH_GPRD
#include "AtomicGPDimer.h"
#endif

class ObjectiveFunction {
protected:
  std::shared_ptr<Matter> matter;
  std::shared_ptr<Parameters> params;

public:
  ObjectiveFunction(std::shared_ptr<Matter> matterPassed,
                    std::shared_ptr<Parameters> paramsPassed)
      : matter{matterPassed}, params{paramsPassed} {}
  virtual ~ObjectiveFunction() {}
  virtual double getEnergy() = 0;
  virtual VectorXd getGradient(bool fdstep = false) = 0;
  virtual void setPositions(VectorXd x) = 0;
  virtual VectorXd getPositions() = 0;
  virtual int degreesOfFreedom() = 0;
  virtual bool isConverged() = 0;
  virtual double getConvergence() = 0;
  virtual VectorXd difference(VectorXd a, VectorXd b) = 0;
};

#endif
