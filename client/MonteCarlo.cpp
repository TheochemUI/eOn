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
#include "MonteCarlo.h"

#include <cmath>
#include <iostream>

using namespace eonc::helpers;

void MonteCarlo::run(int numSteps, double temperature, double stepSize) {

  if (!eonc::io::io_ok(matter->matter2con("movie.con"))) {
    QUILL_LOG_WARNING(log, "Failed to write movie.con header frame");
  }
  AtomMatrix current;
  AtomMatrix trial;
  double T = temperature;
  T = 100.0;
  numSteps = 100;
  stepSize = 0.1;
  int accepts = 0;
  for (int steps = 0; steps < numSteps; steps++) {
    current = matter->getPositions();
    trial = current;
    double ecurrent, etrial;
    ecurrent = matter->getPotentialEnergy();
    if (!eonc::io::io_ok(matter->matter2con("movie.con", true))) {
      QUILL_LOG_WARNING(log, "Failed to append movie.con frame");
    }
    for (int i = 0; i < trial.rows(); i++) {
      for (int j = 0; j < 3; j++) {
        double d = gaussRandom(0.0, stepSize);
        trial(i, j) += d;
      }
    }
    matter->setPositions(trial);
    etrial = matter->getPotentialEnergy();
    double de = ecurrent - etrial;
    QUILL_LOG_INFO(log, "de={}", de);
    if (de <= 0.0) {
      QUILL_LOG_INFO(log, "{}: accept de <= 0.0", steps);
      accepts++;
      continue;
    }
    double r = randomDouble();
    double kB = params.constants.kB;
    double arg = -de / (kB * T);
    QUILL_LOG_DEBUG(log, "arg: {}\n", arg);
    if (arg < -50.0) {
      matter->setPositions(current);
      QUILL_LOG_DEBUG(log, "{}: reject small arg\n", steps);
      continue;
    }

    double p = std::exp(arg);
    if (r < p) {
      QUILL_LOG_DEBUG(log, "{}: accept r<p\n", steps);
      accepts++;
      continue;
    } else {
      matter->setPositions(current);
      QUILL_LOG_DEBUG(log, "{}: reject\n", steps);
    }
  }
  QUILL_LOG_INFO(log, "accepts: {}", accepts);
}
