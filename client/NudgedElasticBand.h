#ifndef NudgedElasticBand_H
#define NudgedElasticBand_H

#include <cmath>
#include <math.h>

#include "Eigen.h"

#include "HelperFunctions.h"
#include "Matter.h"
#include "ObjectiveFunction.h"
#include "Parameters.h"

// NEB method for determining a minimum energy path between two matter objects
class NudgedElasticBand {

public:
  enum class NEBStatus {
    STATUS_GOOD = 0,
    STATUS_INIT = 1,
    STATUS_BAD_MAX_ITERATIONS = 2
  };

  NudgedElasticBand(std::shared_ptr<Matter> initialPassed,
                    std::shared_ptr<Matter> finalPassed,
                    std::shared_ptr<Parameters> parametersPassed,
                    std::shared_ptr<Potential> potPassed);
  ~NudgedElasticBand() = default;

  NudgedElasticBand::NEBStatus compute(void);
  void updateForces(void);
  double convergenceForce(void);
  void findExtrema(void);
  void printImageData(bool writeToFile = false);

  int atoms;
  long numImages, climbingImage, numExtrema;
  std::vector<std::shared_ptr<Matter>> path; // NEB images
  std::vector<std::shared_ptr<AtomMatrix>> tangent;
  std::vector<std::shared_ptr<AtomMatrix>> projectedForce;
  std::vector<double> extremumEnergy;
  std::vector<double> extremumPosition;
  std::vector<double> extremumCurvature;

  long maxEnergyImage;
  bool movedAfterForceCall;

private:
  std::shared_ptr<Parameters> params;
  std::shared_ptr<Potential> pot;
};

class NEBObjectiveFunction : public ObjectiveFunction {
public:
  NEBObjectiveFunction(NudgedElasticBand *nebPassed,
                       std::shared_ptr<Parameters> parametersPassed)
      : ObjectiveFunction(nullptr, parametersPassed), neb{nebPassed} {}
  // This is the odd one out, doesn't take a Matter so we null it

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
};

namespace helper_functions {
namespace neb_paths {
std::vector<Matter> linearPath(const Matter &initImg, const Matter &finalImg,
                               const size_t nimgs);
}
} // namespace helper_functions

#endif
