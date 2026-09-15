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
#include "eon/JobRegistry.h"
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

#include <memory>

namespace eonc {
// Referenced from makeJob so MSVC does not drop this TU from a static
// eonclib (anonymous-namespace registrars would never run).
void forceJobRegistration() {}
} // namespace eonc

#ifdef WITH_GP_SURROGATE
#include "eon/GPSurrogateJob.h"
#endif

namespace {

#define EON_REG_JOB(TYPE, CLS)                                                 \
  const bool eon_reg_##CLS = [] {                                              \
    eonc::registerJob(eonc::JobType::TYPE,                                     \
                      [](std::unique_ptr<eonc::Parameters> p) {                \
                        return std::make_unique<eonc::CLS>(std::move(p));      \
                      });                                                      \
    return true;                                                               \
  }()

EON_REG_JOB(Process_Search, ProcessSearchJob);
EON_REG_JOB(Saddle_Search, SaddleSearchJob);
EON_REG_JOB(Minimization, MinimizationJob);
EON_REG_JOB(Point, PointJob);
EON_REG_JOB(Parallel_Replica, ParallelReplicaJob);
EON_REG_JOB(Safe_Hyperdynamics, SafeHyperJob);
EON_REG_JOB(TAD, TADJob);
EON_REG_JOB(Replica_Exchange, ReplicaExchangeJob);
EON_REG_JOB(Basin_Hopping, BasinHoppingJob);
EON_REG_JOB(Hessian, HessianJob);
EON_REG_JOB(Finite_Difference, FiniteDifferenceJob);
EON_REG_JOB(Nudged_Elastic_Band, NudgedElasticBandJob);
EON_REG_JOB(Dynamics, DynamicsJob);
EON_REG_JOB(Prefactor, PrefactorJob);
EON_REG_JOB(Global_Optimization, GlobalOptimizationJob);
EON_REG_JOB(Structure_Comparison, StructureComparisonJob);
EON_REG_JOB(Monte_Carlo, MonteCarloJob);
#ifdef WITH_GP_SURROGATE
EON_REG_JOB(GP_Surrogate, GPSurrogateJob);
#endif
EON_REG_JOB(OH_TST, OHTSTJob);
EON_REG_JOB(Test, TestJob);

#undef EON_REG_JOB

} // namespace
