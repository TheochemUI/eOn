#include "Parameters.h"
#include "Job.h"
#include "ProcessSearchJob.h"
#include "SaddleSearchJob.h"
#include "MinimizationJob.h"
#include "PointJob.h"
#include "HessianJob.h"
#include "ParallelReplicaJob.h"
#include "SafeHyperJob.h"
#include "TADJob.h"
#include "ReplicaExchangeJob.h"
#include "BasinHoppingJob.h"
#include "FiniteDifferenceJob.h"
#include "NudgedElasticBandJob.h"
#include "DynamicsJob.h"
#include "PrefactorJob.h"
#include "TestJob.h"
#include "GlobalOptimizationJob.h"
#include "StructureComparisonJob.h"
#include "MonteCarloJob.h"

namespace helper_functions {
std::unique_ptr<Job> makeJob(std::unique_ptr<Parameters> params) {
    switch (params->job) {
      case JobType::ProcessSearch: {
          return (std::make_unique<ProcessSearchJob>(std::move(params)));
          break;
      }
      default:
          throw std::runtime_error("No known job could be constructed");
        break;
    }
  }

} // namespace helper_functions
