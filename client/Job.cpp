/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
**
** Copyright (c) 2010--present, eOn Development Team
** All rights reserved.
**
** Repo:
** https://github.com/TheochemUI/eOn
*/
#include "Job.hpp"
// #include "BasinHoppingJob.h"
// #include "DynamicsJob.h"
// #include "FiniteDifferenceJob.h"
// #include "GlobalOptimizationJob.h"
// #include "HessianJob.h"
// #include "MinimizationJob.h"
// #include "MonteCarloJob.h"
// #include "NudgedElasticBandJob.h"
// #include "ParallelReplicaJob.h"
#include "BaseStructures.h"
#include "MinimizationJob.hpp"
#include "Optimizer.h"
#include "Parser.hpp"
#include "PointJob.h"
#include "matter/StructComparer.hpp"
#include "parsers/ParseJob.hpp"
#include "parsers/ParseOptim.hpp"
// #include "PrefactorJob.h"
// #include "ProcessSearchJob.h"
// #include "ReplicaExchangeJob.h"
// #include "SaddleSearchJob.h"
// #include "SafeHyperJob.h"
// #include "StructureComparisonJob.h"
// #include "TADJob.h"
// #include "TestJob.h"

#ifdef WITH_GP_SURROGATE
#include "GPSurrogateJob.h"
#endif

namespace eonc {

JobVariant mkJob(const toml::table &config) {
  config_section(config, "Main");
  auto jtype = get_enum_toml<JobType>(config["Main"]["job"]).value();
  switch (jtype) {
  case JobType::Point: {
    return PointJob();
  }
  case JobType::Minimization: {
    config_section(config, "Optimizer");
    auto params = eonc::OptimBase::Params();
    eonc::opt::from_toml(params, config["Optimizer"]);
    return MinimizationJob(params);
  }
  case JobType::Structure_Comparison: {
    config_section(config, "Structure_Comparison");
    auto params = mat::StructComparer::Params();
    eonc::job::from_toml(params, config["Structure_Comparison"]);
    auto sc = eonc::mat::StructComparer(params);
    return StructureComparisonJob(sc);
  }
  default: {
    throw std::runtime_error("No known job could be constructed from mkJob");
  }
  }
}

bool runJob(JobVariant &job, std::vector<Matter> &mats) {
  bool result = false;
  if (std::holds_alternative<StructureComparisonJob>(job) && mats.size() >= 2) {
    result = eonc::JobRunner(job, mats[0], mats[1]);
  } else if (std::holds_alternative<PointJob>(job) && !mats.empty()) {
    result = eonc::JobRunner(job, mats[0]);
  } else if (std::holds_alternative<MinimizationJob>(job) && !mats.empty()) {
    result = eonc::JobRunner(job, mats[0]);
  } else {
    throw std::runtime_error("Invalid job type or arguments");
  }
  return result;
}

// std::unique_ptr<JobBase>
// makeJob(const toml::table &config,
//         std::optional<std::reference_wrapper<Matter>> mat) {
//   config_section(config, "Main");
//   auto jtype = get_enum_toml<JobType>(config["Main"]["job"]);
//   switch (jtype) {
//   // case JobType::Process_Search: {
//   //   return (std::make_unique<ProcessSearchJob>(std::move(params)));
//   //   break;
//   // }
//   // case JobType::Saddle_Search: {
//   //   return (std::make_unique<SaddleSearchJob>(std::move(params)));
//   //   break;
//   // }
//   // case JobType::Minimization: {
//   //   return (std::make_unique<MinimizationJob>(std::move(params)));
//   //   break;
//   // }
//   case JobType::Point: {
//     return (std::make_unique<PointJob>(mat->get()));
//     break;
//   }
//     //   case JobType::Parallel_Replica: {
//     //     return (std::make_unique<ParallelReplicaJob>(std::move(params)));
//     //     break;
//     //   }
//     //   case JobType::Safe_Hyperdynamics: {
//     //     return (std::make_unique<SafeHyperJob>(std::move(params)));
//     //     break;
//     //   }
//     //   case JobType::TAD: {
//     //     return (std::make_unique<TADJob>(std::move(params)));
//     //     break;
//     //   }
//     //   case JobType::Replica_Exchange: {
//     //     return (std::make_unique<ReplicaExchangeJob>(std::move(params)));
//     //     break;
//     //   }
//     //   case JobType::Basin_Hopping: {
//     //     return (std::make_unique<BasinHoppingJob>(std::move(params)));
//     //     break;
//     //   }
//     //   case JobType::Hessian: {
//     //     return (std::make_unique<HessianJob>(std::move(params)));
//     //     break;
//     //   }
//     //   case JobType::Finite_Difference: {
//     //     return (std::make_unique<FiniteDifferenceJob>(std::move(params)));
//     //     break;
//     //   }
//     //   case JobType::Nudged_Elastic_Band: {
//     //     return
//     (std::make_unique<NudgedElasticBandJob>(std::move(params)));
//     //     break;
//     //   }
//     //   case JobType::Dynamics: {
//     //     return (std::make_unique<DynamicsJob>(std::move(params)));
//     //     break;
//     //   }
//     //   case JobType::Prefactor: {
//     //     return (std::make_unique<PrefactorJob>(std::move(params)));
//     //     break;
//     //   }
//     //   case JobType::Global_Optimization: {
//     //     return
//     (std::make_unique<GlobalOptimizationJob>(std::move(params)));
//     //     break;
//     //   }
//     //   case JobType::Structure_Comparison: {
//     //     return
//     (std::make_unique<StructureComparisonJob>(std::move(params)));
//     //     break;
//     //   }
//     //   case JobType::Monte_Carlo: {
//     //     return (std::make_unique<MonteCarloJob>(std::move(params)));
//     //     break;
//     //   }
//     // #ifdef WITH_GP_SURROGATE
//     //   case JobType::GP_Surrogate: {
//     //     return (std::make_unique<GPSurrogateJob>(std::move(params)));
//     //     break;
//     //   }
//     // #endif
//   default:
//     throw std::runtime_error("No known job could be constructed");
//     break;
//   }
// }

} // namespace eonc
