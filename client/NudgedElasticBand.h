#ifndef NudgedElasticBand_H
#define NudgedElasticBand_H

#include <cmath>
#include <filesystem>
#include <math.h>

#include "Eigen.h"

#include "HelperFunctions.h"
#include "LowestEigenmode.h"
#include "Matter.h"
#include "ObjectiveFunction.h"
#include "Parameters.h"

// NEB method for determining a minimum energy path between two matter objects
class NudgedElasticBand {

public:
  enum class NEBStatus {
    GOOD = 0,
    INIT = 1,
    BAD_MAX_ITERATIONS = 2,
    RUNNING,
    MAX_UNCERTAINITY
  };
  NudgedElasticBand(std::shared_ptr<Matter> initialPassed,
                    std::shared_ptr<Matter> finalPassed,
                    std::shared_ptr<Parameters> parametersPassed,
                    std::shared_ptr<Potential> potPassed);
  NudgedElasticBand(std::vector<Matter> initPath,
                    std::shared_ptr<Parameters> parametersPassed,
                    std::shared_ptr<Potential> potPassed);
  ~NudgedElasticBand() = default;

  NudgedElasticBand::NEBStatus compute(void);
  NudgedElasticBand::NEBStatus getStatus() { return this->status; };
  void updateForces(void);
  double convergenceForce(void);
  void findExtrema(void);
  void printImageData(bool writeToFile = false, size_t idx = 0);
  std::vector<std::shared_ptr<LowestEigenmode>> eigenmode_solvers;

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
  double ksp;
  double k_u;   // Upper-bound value for the spring constant
  double k_l;   // Lower-bound value for the spring constant
  double E_ref; // Reference energy chosen to be equal to the max energy of the
                // reactant or product energy minimum

private:
  std::shared_ptr<Parameters> params;
  std::shared_ptr<Potential> pot;
  NEBStatus status;
  std::shared_ptr<spdlog::logger> log;
};

class NEBObjectiveFunction : public ObjectiveFunction {
public:
  NEBObjectiveFunction(NudgedElasticBand *nebPassed,
                       std::shared_ptr<Parameters> parametersPassed)
      : ObjectiveFunction(nullptr, parametersPassed),
        neb{nebPassed} {}
  // This is the odd one out, doesn't take a Matter so we null it

  ~NEBObjectiveFunction(void) {};

  VectorXd getGradient(bool fdstep = false);
  double getEnergy();
  void setPositions(VectorXd x);
  VectorXd getPositions();
  int degreesOfFreedom();
  bool isConverged();
  bool isUncertain();
  double getConvergence();
  VectorXd difference(VectorXd a, VectorXd b);
  NudgedElasticBand::NEBStatus status;

private:
  NudgedElasticBand *neb;
};

namespace helper_functions {
namespace neb_paths {
std::vector<Matter> linearPath(const Matter &initImg, const Matter &finalImg,
                               const size_t nimgs);
std::vector<Matter>
filePathInit(const std::vector<std::filesystem::path> &fsrcs,
             const Matter &refImg, const size_t nimgs);
/**
 * @brief Reads a file where each line contains a path to another file.
 *
 * @param listFilePath The path to the file containing the list of file paths.
 * @return A vector of filesystem paths. Returns an empty vector if the
 * file cannot be opened.
 */
std::vector<std::filesystem::path>
readFilePaths(const std::string &listFilePath);
} // namespace neb_paths
} // namespace helper_functions

#endif
