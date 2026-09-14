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
#include "eon/Job.h"
#include "eon/BasinHoppingJob.h"
#include "eon/DynamicsJob.h"
#include "eon/FiniteDifferenceJob.h"
#include "eon/GlobalOptimizationJob.h"
#include "eon/HessianJob.h"
#include "eon/MinimizationJob.h"
#include "eon/MonteCarloJob.h"
#include "eon/NudgedElasticBandJob.h"
#include "eon/OHTSTJob.h"
#include "eon/ParallelReplicaJob.h"
#include "eon/Parameters.h"
#include "eon/PointJob.h"
#include "eon/PrefactorJob.h"
#include "eon/ProcessSearchJob.h"
#include "eon/ReplicaExchangeJob.h"
#include "eon/SaddleSearchJob.h"
#include "eon/SafeHyperJob.h"
#include "eon/StructureComparisonJob.h"
#include "eon/TADJob.h"
#include "eon/TestJob.h"

#ifdef WITH_GP_SURROGATE
#include "eon/GPSurrogateJob.h"
#endif

namespace eonc::helpers {
std::unique_ptr<Job> makeJob(std::unique_ptr<Parameters> params) {
  switch (params->main_options.job) {
    using enum JobType;
  case Process_Search: {
    return (std::make_unique<ProcessSearchJob>(std::move(params)));
    break;
  }
  case Saddle_Search: {
    return (std::make_unique<SaddleSearchJob>(std::move(params)));
    break;
  }
  case Minimization: {
    return (std::make_unique<MinimizationJob>(std::move(params)));
    break;
  }
  case Point: {
    return (std::make_unique<PointJob>(std::move(params)));
    break;
  }
  case Parallel_Replica: {
    return (std::make_unique<ParallelReplicaJob>(std::move(params)));
    break;
  }
  case Safe_Hyperdynamics: {
    return (std::make_unique<SafeHyperJob>(std::move(params)));
    break;
  }
  case TAD: {
    return (std::make_unique<TADJob>(std::move(params)));
    break;
  }
  case Replica_Exchange: {
    return (std::make_unique<ReplicaExchangeJob>(std::move(params)));
    break;
  }
  case Basin_Hopping: {
    return (std::make_unique<BasinHoppingJob>(std::move(params)));
    break;
  }
  case Hessian: {
    return (std::make_unique<HessianJob>(std::move(params)));
    break;
  }
  case Finite_Difference: {
    return (std::make_unique<FiniteDifferenceJob>(std::move(params)));
    break;
  }
  case Nudged_Elastic_Band: {
    return (std::make_unique<NudgedElasticBandJob>(std::move(params)));
    break;
  }
  case Dynamics: {
    return (std::make_unique<DynamicsJob>(std::move(params)));
    break;
  }
  case Prefactor: {
    return (std::make_unique<PrefactorJob>(std::move(params)));
    break;
  }
  case Global_Optimization: {
    return (std::make_unique<GlobalOptimizationJob>(std::move(params)));
    break;
  }
  case Structure_Comparison: {
    return (std::make_unique<StructureComparisonJob>(std::move(params)));
    break;
  }
  case Monte_Carlo: {
    return (std::make_unique<MonteCarloJob>(std::move(params)));
    break;
  }
#ifdef WITH_GP_SURROGATE
  case GP_Surrogate: {
    return (std::make_unique<GPSurrogateJob>(std::move(params)));
    break;
  }
#endif
  case OH_TST: {
    return (std::make_unique<OHTSTJob>(std::move(params)));
  }
  case Test: {
    return (std::make_unique<TestJob>(std::move(params)));
  }
  default:
    throw std::runtime_error("No known job could be constructed");
    break;
  }
}

} // namespace eonc::helpers
