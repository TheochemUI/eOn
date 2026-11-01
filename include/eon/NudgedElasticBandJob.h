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
#pragma once
#include "EonLogger.h"

#include "Job.h"
#include "Matter.h"
#include "NudgedElasticBand.h"
#include "Parameters.h"
#include "PrfoRefine.h"

#include <optional>

namespace eonc {

class NudgedElasticBandJob : public Job {

public:
  NudgedElasticBandJob(std::unique_ptr<Parameters> parameters, Runtime &rt)
      : Job(std::move(parameters), rt), fCallsNEB{0} {}
  ~NudgedElasticBandJob(void) = default;
  std::vector<std::string> run(void);

private:
  // functions
  void printEndState(NudgedElasticBand::NEBStatus status);
  void saveData(NudgedElasticBand::NEBStatus status, NudgedElasticBand *neb);
  void refineClimbingImage(NudgedElasticBand &neb);

  // variables
  std::vector<std::string> returnFiles;
  size_t fCallsNEB;
  std::optional<PrfoRefineResult> refined_;
  size_t fCallsRefine_{0};
  eonc::log::Scoped m_log;
};

} // namespace eonc
