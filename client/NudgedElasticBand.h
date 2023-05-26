#ifndef NudgedElasticBand_H
#define NudgedElasticBand_H

#include <cmath>
#include <math.h>

#include "Eigen.h"

#include "HelperFunctions.h"
#include "Matter.h"
#include "ObjectiveFunction.h"
#include "Parameters.h"

class Matter;
class Parameters;

// NEB method for determining a minimum energy path between two matter objects
class NudgedElasticBand {

public:
  enum class NEBStatus {
    STATUS_GOOD = 0,
    STATUS_INIT = 1,
    STATUS_BAD_MAX_ITERATIONS = 2
  };

  NudgedElasticBand(Matter *initialPassed, Matter *finalPassed,
                    Parameters *parametersPassed);
  ~NudgedElasticBand();

  void clean(void);
  NudgedElasticBand::NEBStatus compute(void);
  void updateForces(void);
  double convergenceForce(void);
  void findExtrema(void);
  void printImageData(bool writeToFile = false);

  int atoms;
  long images, climbingImage, numExtrema;
  Matter **image; // NEB images
  AtomMatrix **tangent;
  AtomMatrix **projectedForce;
  bool movedAfterForceCall;
  double *extremumEnergy;
  double *extremumPosition;
  double *extremumCurvature;

  long maxEnergyImage;

private:
  Parameters *parameters;
};

class NEBObjectiveFunction : public ObjectiveFunction {
public:
  NEBObjectiveFunction(NudgedElasticBand *nebPassed,
                       Parameters *parametersPassed)
      : neb{nebPassed}, parameters{parametersPassed} {}

  ~NEBObjectiveFunction(void){};

  VectorXd getGradient(bool fdstep = false);
  double getEnergy();
  void setPositions(VectorXd x);
  VectorXd getPositions();
  int degreesOfFreedom();
  bool isConverged();
  double getConvergence();
  VectorXd difference(VectorXd a, VectorXd b);

private:
  NudgedElasticBand *neb;
  Parameters *parameters;
};

#endif
