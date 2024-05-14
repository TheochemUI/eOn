#ifndef POTENTIAL_H
#define POTENTIAL_H

#include "Eigen.h"
#include "Parameters.h"
#include <algorithm>
#include <limits>
#include <memory>
#include <optional>

class Potential {
protected:
  PotType ptype;
  std::shared_ptr<Parameters> m_params;
  std::shared_ptr<spdlog::logger> m_log;

public:
  size_t forceCallCounter;
  Potential(PotType a_ptype, std::shared_ptr<Parameters> a_params)
      : ptype{a_ptype},
        forceCallCounter{0},
        m_params{a_params} {
    initializeLogger();
  }
  Potential(std::shared_ptr<Parameters> a_params)
      : ptype{a_params->potential},
        m_params{a_params} {
    initializeLogger();
  }

  virtual ~Potential() {
    std::string potentialName = helper_functions::getPotentialName(getType());
    m_log->info("[{}] called potential {} times", potentialName,
                forceCallCounter);
  }
  static int fcalls;
  static int fcallsTotal;
  static int wu_fcallsTotal;
  static double totalUserTime;

  // Does not take into account the fixed / free atoms
  // Variance here is null when not needed and that's OK
  void virtual force(long nAtoms, const double *positions, const int *atomicNrs,
                     double *forces, double *energy, double *variance,
                     const double *box) = 0;
  std::tuple<double, AtomMatrix>
  get_ef(const AtomMatrix pos, const VectorXi atmnrs, const Matrix3d box);
  PotType getType() { return this->ptype; };
  void initializeLogger() {
    if (!spdlog::get("_potcalls")) {
      // Create logger if it doesn't exist
      m_log = spdlog::basic_logger_mt("_potcalls", "_potcalls.log", true);
      m_log->set_pattern("[%l] [%Y-%m-%d %H:%M:%S] %v");
    } else {
      // Use existing logger
      m_log = spdlog::get("_potcalls");
    }
    std::string potentialName = helper_functions::getPotentialName(getType());
    m_log->info("[{}] created", potentialName);
  }
};

namespace helper_functions {
std::shared_ptr<Potential> makePotential(std::shared_ptr<Parameters> params);
std::shared_ptr<Potential> makePotential(PotType ptype,
                                         std::shared_ptr<Parameters> params);
} // namespace helper_functions

// int Potential::fcalls = 0;
// int Potential::fcallsTotal = 0;
// int Potential::wu_fcallsTotal = 0;
// double Potential::totalUserTime = 0.0;
#endif
