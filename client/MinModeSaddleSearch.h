#ifndef MIN_MODE_SADDLE_SEARCH_H
#define MIN_MODE_SADDLE_SEARCH_H

#include "Eigen.h"
#include "LowestEigenmode.h"
#include "Matter.h"
#include "Optimizer.h"
#include "SaddleSearchMethod.h"

#include <string>

class MinModeSaddleSearch : public SaddleSearchMethod {

public:
  enum {
    // DO NOT CHANGE THE ORDER OF THIS LIST
    STATUS_GOOD,                           // 0
    STATUS_INIT,                           // 1
    STATUS_BAD_NO_CONVEX,                  // 2
    STATUS_BAD_HIGH_ENERGY,                // 3
    STATUS_BAD_MAX_CONCAVE_ITERATIONS,     // 4
    STATUS_BAD_MAX_ITERATIONS,             // 5
    STATUS_BAD_NOT_CONNECTED,              // 6
    STATUS_BAD_PREFACTOR,                  // 7
    STATUS_BAD_HIGH_BARRIER,               // 8
    STATUS_BAD_MINIMA,                     // 9
    STATUS_FAILED_PREFACTOR,               // 10
    STATUS_POTENTIAL_FAILED,               // 11
    STATUS_NONNEGATIVE_ABORT,              // 12
    STATUS_NONLOCAL_ABORT,                 // 13
    STATUS_NEGATIVE_BARRIER,               // 14
    STATUS_BAD_MD_TRAJECTORY_TOO_SHORT,    // 15
    STATUS_BAD_NO_NEGATIVE_MODE_AT_SADDLE, // 16
    STATUS_BAD_NO_BARRIER,                 // 17
    STATUS_ZEROMODE_ABORT,                 // 18
    STATUS_OPTIMIZER_ERROR                 // 19
  };

  MinModeSaddleSearch(std::shared_ptr<Matter> matterPassed,
                      AtomMatrix modePassed, double reactantEnergyPassed,
                      std::shared_ptr<Parameters> parametersPassed,
                      std::shared_ptr<Potential> potPassed);
  ~MinModeSaddleSearch() = default;
  AtomMatrix getEigenvector(); // lowest eigenmode
  double getEigenvalue();      // estimate for the lowest eigenvalue

  int run();

  int forcecalls;
  int iteration;
  int status;

private:
  AtomMatrix mode;
  std::shared_ptr<Matter> matter;
  std::shared_ptr<LowestEigenmode>
      minModeMethod; // shared with the objective func
  double reactantEnergy;
  shared_ptr<spdlog::logger> log;
};

#endif
