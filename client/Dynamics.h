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

#include "HelperFunctions.h"
#include "Matter.h"
#include "Parameters.h"

#include "Eigen.h"
#include "EonLogger.h"

namespace eonc {

/// Configuration extracted from Parameters for Dynamics.
struct DynamicsConfig {
  double time_step{0.0};
  long steps{0};
  std::string thermostat_kind{"none"};
  double andersen_alpha{1.0};
  double andersen_tcol{0.0};
  double nose_mass{1.0};
  double langevin_friction{0.0};
  double kB{8.6173324e-5};
  double timeUnit{10.1805055};
  double temperature{300.0};
  bool write_movies{false};
  long write_movies_interval{1};

  static DynamicsConfig fromParams(const Parameters &p) {
    return {p.dynamics_options.time_step,
            p.dynamics_options.steps,
            p.thermostat_options.kind,
            p.thermostat_options.andersen_alpha,
            p.thermostat_options.andersen_tcol,
            p.thermostat_options.nose_mass,
            p.thermostat_options.langevin_friction,
            p.constants.kB,
            p.constants.timeUnit,
            p.main_options.temperature,
            p.debug_options.write_movies,
            p.debug_options.write_movies_interval};
  }
};

class Dynamics {

public:
  static const char ANDERSEN[];
  static const char NOSE_HOOVER[];
  static const char LANGEVIN[];
  static const char NONE[];

  Dynamics(Matter *matter, const DynamicsConfig &config);

  [[deprecated("Pass DynamicsConfig directly")]]
  Dynamics(Matter *matter, const Parameters &parameters)
      : Dynamics(matter, DynamicsConfig::fromParams(parameters)) {}

  ~Dynamics();

  void setTemperature(double temperature);
  void oneStep(int stepNumber = -1);
  void velocityVerlet();
  void run();
  void andersenCollision();
  void setThermalVelocity();
  void rescaleVelocity();
  void noseHooverVerlet();
  void langevinVerlet();

private:
  long nAtoms{0}, nFreeCoords{0};

  Matter *matter;
  DynamicsConfig m_config;

  double dt{0.0};
  double kB{0.0};
  double temperature{0.0};
  double vxi1{0.0}, vxi2{0.0}, xi1{0.0}, xi2{0.0};
  eonc::log::Scoped log;
};

} // namespace eonc

using eonc::Dynamics;
using eonc::DynamicsConfig;
