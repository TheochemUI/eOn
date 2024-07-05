#pragma once

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
      : matter{matterPassed},
        params{paramsPassed} {}
  virtual ~ObjectiveFunction() {}
  virtual double getEnergy() = 0;
  virtual VectorType getGradient(bool fdstep = false) = 0;
  virtual void setPositions(VectorType x) = 0;
  virtual VectorType getPositions() = 0;
  virtual int degreesOfFreedom() = 0;
  virtual bool isConverged() = 0;
  virtual double getConvergence() = 0;
  virtual VectorType difference(VectorType a, VectorType b) = 0;
};
