#include "BaseStructures.h"
#include <cassert>
#include <iostream>

namespace helper_functions {

std::string getJobName(JobType jtype) {
  switch (jtype) {
  case JobType::ProcessSearch: {
    return "process_search"s;
    break;
  }
  case JobType::SaddleSearch: {
    return "saddle_search"s;
    break;
  }
  case JobType::Minimization: {
    return "minimization"s;
    break;
  }
  case JobType::Point: {
    return "point"s;
    break;
  }
  case JobType::ParallelReplica: {
    return "parallel_replica"s;
    break;
  }
  case JobType::SafeHyperdynamics: {
    return "safe_hyperdynamics"s;
    break;
  }
  case JobType::TAD: {
    return "tad"s;
    break;
  }
  case JobType::ReplicaExchange: {
    return "replica_exchange"s;
    break;
  }
  case JobType::BasinHopping: {
    return "basin_hopping"s;
    break;
  }
  case JobType::Hessian: {
    return "hessian"s;
    break;
  }
  case JobType::FiniteDifference: {
    return "finite_difference"s;
    break;
  }
  case JobType::NEB: {
    return "nudged_elastic_band"s;
    break;
  }
  case JobType::Dynamics: {
    return "dynamics"s;
    break;
  }
  case JobType::Prefactor: {
    return "prefactor"s;
    break;
  }
  case JobType::GlobalOptimization: {
    return "global_optimization"s;
    break;
  }
  case JobType::StructureComparison: {
    return "structure_comparison"s;
    break;
  }
  case JobType::MonteCarlo: {
    return "monte_carlo"s;
    break;
  }
  case JobType::GPSurrogate: {
    return "gp_surrogate"s;
    break;
  }
  default:
    return "unknown job"s;
    break;
  }
}

JobType getJobType(const std::string jname) {
  if (jname == "process_search"s) {
    return JobType::ProcessSearch;
  } else if (jname == "saddle_search"s) {
    return JobType::SaddleSearch;
  } else if (jname == "minimization"s) {
    return JobType::Minimization;
  } else if (jname == "point"s) {
    return JobType::Point;
  } else if (jname == "parallel_replica"s) {
    return JobType::ParallelReplica;
  } else if (jname == "safe_hyperdynamics"s) {
    return JobType::SafeHyperdynamics;
  } else if (jname == "tad"s) {
    return JobType::TAD;
  } else if (jname == "replica_exchange"s) {
    return JobType::ReplicaExchange;
  } else if (jname == "basin_hopping"s) {
    return JobType::BasinHopping;
  } else if (jname == "hessian"s) {
    return JobType::Hessian;
  } else if (jname == "finite_difference"s) {
    return JobType::FiniteDifference;
  } else if (jname == "nudged_elastic_band"s) {
    return JobType::NEB;
  } else if (jname == "dynamics"s) {
    return JobType::Dynamics;
  } else if (jname == "prefactor"s) {
    return JobType::Prefactor;
  } else if (jname == "global_optimization"s) {
    return JobType::GlobalOptimization;
  } else if (jname == "structure_comparison"s) {
    return JobType::StructureComparison;
  } else if (jname == "monte_carlo"s) {
    return JobType::MonteCarlo;
  } else if (jname == "gp_surrogate"s) {
    return JobType::GPSurrogate;
  } else {
    return JobType::Unknown;
  }
}

std::string getRunStatusName(RunStatus rstype) {
  switch (rstype) {
  case RunStatus::GOOD: {
    return "PASS: good"s;
    break;
  }
  case RunStatus::MAX_ITERATIONS: {
    return "FAIL: max iterations reached"s;
    break;
  }
  case RunStatus::POTENTIAL_FAILED: {
    return "FAIL: potential failed"s;
    break;
  }
  default:
    return "unknown"s;
    break;
  }
}

std::string getOptName(OptType a_otype) {
  switch (a_otype) {
  case OptType::None: {
    return "No optimizer"s;
  }
  case OptType::Unknown: {
    return "Unknown optimizer"s;
  }
  case OptType::QuickMin: {
    return "QuickMin optimizer"s;
  }
  case OptType::ConjugateGradient: {
    return "Conjugate Gradient optimizer"s;
  }
  case OptType::LBFGS: {
    return "Limited memory Broyden-Fletcher-Goldfarb-Shanno optimizer"s;
  }
  case OptType::FIRE: {
    return "Fast Inertial Relaxation Engine optimizer"s;
  }
  case OptType::SteepestDescent: {
    return "Steepest Descent optimizer"s;
  }
  default:
    throw std::runtime_error("[Error] Invalid Optimizer type!!");
  }
}
OptType getOptType(std::string a_oname) {
  if (a_oname == "cg"s) {
    return OptType::ConjugateGradient;
  } else if (a_oname == "qm"s) {
    return OptType::QuickMin;
  } else if (a_oname == "lbfgs"s) {
    return OptType::LBFGS;
  } else if (a_oname == "fire"s) {
    return OptType::FIRE;
  } else if (a_oname == "sd"s) {
    return OptType::SteepestDescent;
  } else if (a_oname == "none"s) {
    return OptType::None;
  } else if (a_oname == "unknown"s) {
    return OptType::Unknown;
  } else {
    throw std::runtime_error("[Error] Invalid Optimizer string!!");
  }
}
} // namespace helper_functions
