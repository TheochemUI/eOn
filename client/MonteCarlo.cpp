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
#include "eon/MonteCarlo.h"

#include <cmath>
#include <stdexcept>

namespace eonc {

void MonteCarlo::run(int numSteps, double temperature, double stepSize) {
  if (numSteps <= 0) {
    throw std::invalid_argument("MonteCarlo: steps must be positive");
  }
  if (!(temperature > 0.0)) {
    throw std::invalid_argument("MonteCarlo: temperature must be positive");
  }
  if (!(stepSize > 0.0)) {
    throw std::invalid_argument("MonteCarlo: step_size must be positive");
  }

  if (!eonc::io::io_ok(matter->matter2con("movie.con"))) {
    QUILL_LOG_WARNING(log, "Failed to write movie.con header frame");
  }

  const double kB = params.constants().kB;
  int accepts = 0;
  for (int steps = 0; steps < numSteps; ++steps) {
    const AtomMatrix current = matter->getPositions();
    const double ecurrent = matter->getPotentialEnergy();
    if (!eonc::io::io_ok(matter->matter2con("movie.con", true))) {
      QUILL_LOG_WARNING(log, "Failed to append movie.con frame");
    }

    AtomMatrix trial = current;
    for (int i = 0; i < trial.rows(); ++i) {
      if (matter->getFixed(i)) {
        continue;
      }
      for (int j = 0; j < 3; ++j) {
        trial(i, j) += eonc::rng::gaussRandom(0.0, stepSize);
      }
    }
    matter->setPositions(trial);
    const double etrial = matter->getPotentialEnergy();
    const double de = etrial - ecurrent;
    bool accept = de <= 0.0;
    if (!accept) {
      const double arg = -de / (kB * temperature);
      accept = (arg >= -50.0) && (eonc::rng::randomDouble() < std::exp(arg));
    }
    if (accept) {
      ++accepts;
    } else {
      matter->setPositions(current);
    }
  }
  QUILL_LOG_INFO(log, "accepts: {}", accepts);
}

} // namespace eonc
