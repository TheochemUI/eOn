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

/// Fixed-energy line-integral nudged elastic band (LI-NEB) through
/// rgsaddle: Ásgeirsson, Arnaldsson and Jónsson, J. Chem. Phys. 148,
/// 102334 (2018). Needs `-Dwith_rgsaddle=true` against a rgsaddle that
/// exports the lineb C API (`rgsaddle_lineb_solve`). The built-in
/// ring-polymer Instanton job is unchanged.
class LinebJob : public Job {
public:
  LinebJob(std::unique_ptr<Parameters> parameters, Runtime &rt)
      : Job(std::move(parameters), rt) {}
  LinebJob(std::shared_ptr<Potential> potPassed, const Parameters &parameters)
      : Job(std::move(potPassed), parameters) {}
  ~LinebJob(void) = default;
  std::vector<std::string> run(void) override;
};

/// True when this build linked a rgsaddle that exports the lineb C API.
bool linebAvailable();

} // namespace eonc
