#include "Job.h"
#include "BasinHoppingJob.h"
#include "DynamicsJob.h"
#include "FiniteDifferenceJob.h"
#include "GlobalOptimizationJob.h"
#include "HessianJob.h"
#include "MinimizationJob.h"
#include "MonteCarloJob.h"
#include "NudgedElasticBandJob.h"
#include "ParallelReplicaJob.h"
#include "Parameters.h"
#include "PointJob.h"
#include "PrefactorJob.h"
#include "ProcessSearchJob.h"
#include "ReplicaExchangeJob.h"
#include "SaddleSearchJob.h"
#include "SafeHyperJob.h"
#include "StructureComparisonJob.h"
#include "TADJob.h"
#include "TestJob.h"

#ifdef WITH_GP_SURROGATE
#include "GPSurrogateJob.h"
#endif

namespace helper_functions {
std::unique_ptr<Job> makeJob(std::unique_ptr<Parameters> params) {
  switch (params->job) {
  case JobType::ProcessSearch: {
    return (std::make_unique<ProcessSearchJob>(std::move(params)));
    break;
  }
  case JobType::SaddleSearch: {
    return (std::make_unique<SaddleSearchJob>(std::move(params)));
    break;
  }
  case JobType::Minimization: {
    return (std::make_unique<MinimizationJob>(std::move(params)));
    break;
  }
  case JobType::Point: {
    return (std::make_unique<PointJob>(std::move(params)));
    break;
  }
  case JobType::ParallelReplica: {
    return (std::make_unique<ParallelReplicaJob>(std::move(params)));
    break;
  }
  case JobType::SafeHyperdynamics: {
    return (std::make_unique<SafeHyperJob>(std::move(params)));
    break;
  }
  case JobType::TAD: {
    return (std::make_unique<TADJob>(std::move(params)));
    break;
  }
  case JobType::ReplicaExchange: {
    return (std::make_unique<ReplicaExchangeJob>(std::move(params)));
    break;
  }
  case JobType::BasinHopping: {
    return (std::make_unique<BasinHoppingJob>(std::move(params)));
    break;
  }
  case JobType::Hessian: {
    return (std::make_unique<HessianJob>(std::move(params)));
    break;
  }
  case JobType::FiniteDifference: {
    return (std::make_unique<FiniteDifferenceJob>(std::move(params)));
    break;
  }
  case JobType::NEB: {
    return (std::make_unique<NudgedElasticBandJob>(std::move(params)));
    break;
  }
  case JobType::Dynamics: {
    return (std::make_unique<DynamicsJob>(std::move(params)));
    break;
  }
  case JobType::Prefactor: {
    return (std::make_unique<PrefactorJob>(std::move(params)));
    break;
  }
  case JobType::GlobalOptimization: {
    return (std::make_unique<GlobalOptimizationJob>(std::move(params)));
    break;
  }
  case JobType::StructureComparison: {
    return (std::make_unique<StructureComparisonJob>(std::move(params)));
    break;
  }
  case JobType::MonteCarlo: {
    return (std::make_unique<MonteCarloJob>(std::move(params)));
    break;
  }
#ifdef WITH_GP_SURROGATE
  case JobType::GPSurrogate: {
    return (std::make_unique<GPSurrogateJob>(std::move(params)));
    break;
  }
#endif
  default:
    throw std::runtime_error("No known job could be constructed");
    break;
  }
}

} // namespace helper_functions
