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
#include "Job.h"
#include "Parameters.h"

namespace eonc {

/// Ring-polymer instanton ([Instanton]): the tunnelling splitting between
/// two minima beyond one-dimensional WKB (mode splitting), or the thermal
/// rate through a saddle (mode rate): the ring below the crossover
/// temperature, and the parabolic barrier factor times harmonic
/// transition-state theory above it, written to instanton.con and
/// results.dat.
class InstantonJob : public Job {
public:
  InstantonJob(std::unique_ptr<Parameters> parameters, Runtime &rt)
      : Job(std::move(parameters), rt) {}
  ~InstantonJob(void) = default;
  std::vector<std::string> run(void) override;
};

/// External orbit. `evaluate` returns 0 and writes `action` when it owns
/// the rate. A null `evaluate` leaves the built-in ring polymer in place.
struct ExternalInstanton {
  int (*evaluate)(const double *position, int n_coord, const double *masses,
                  int n_atom, double temperature_k, double *action,
                  void *user);
  void *user;
};

void set_external_instanton(ExternalInstanton hook);
ExternalInstanton external_instanton();

} // namespace eonc
