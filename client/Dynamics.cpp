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
#include "Dynamics.h"
#include "EonLogger.h"

#include <cmath>

using namespace eonc::helpers;

const char Dynamics::ANDERSEN[] = "andersen";
const char Dynamics::NOSE_HOOVER[] = "nose_hoover";
const char Dynamics::LANGEVIN[] = "langevin";
const char Dynamics::NONE[] = "none";

Dynamics::Dynamics(Matter *matter_in, const DynamicsConfig &config)
    : matter{matter_in},
      m_config{config} {
  dt = m_config.time_step;
  nAtoms = matter->numberOfAtoms();
  nFreeCoords = matter->numberOfFreeAtoms() * 3;
  temperature = m_config.temperature;
  kB = m_config.kB;
  vxi1 = vxi2 = xi1 = xi2 = 0.0;
}

Dynamics::~Dynamics() = default;

void Dynamics::setTemperature(double temperature_in) {
  temperature = temperature_in;
}

void Dynamics::oneStep(int stepNumber) {
  if (m_config.thermostat_kind == ANDERSEN) {
    andersenCollision();
    velocityVerlet();
  } else if (m_config.thermostat_kind == NOSE_HOOVER) {
    noseHooverVerlet();
  } else if (m_config.thermostat_kind == LANGEVIN) {
    langevinVerlet();
  } else if (m_config.thermostat_kind == NONE) {
    velocityVerlet();
  }

  if (stepNumber != -1) {
    if (stepNumber == 1) {
      QUILL_LOG_DEBUG(log, "{} {:8s} {:10s} {:12s} {:12s} {:10s}\n",
                      "[Dynamics]", "Step", "KE", "PE", "TE", "KinT");
    }
    double kinE = matter->getKineticEnergy();
    double potE = matter->getPotentialEnergy();
    double kinT = 2.0 * kinE / nFreeCoords / kB;

    if (stepNumber % m_config.write_movies_interval == 0) {
      QUILL_LOG_DEBUG(log, "{} {:8} {:10.4} {:12.4} {:12.4} {:10.2}\n",
                      "[Dynamics]", stepNumber, kinE, potE, kinE + potE, kinT);
    }
  }
}

void Dynamics::velocityVerlet() {
  AtomMatrix positions = matter->getPositions();
  AtomMatrix velocities = matter->getVelocities();
  AtomMatrix accInit = matter->getAccelerations();

  positions += dt * velocities + 0.5 * dt * dt * accInit;
  matter->setPositions(positions);

  AtomMatrix accFinal = matter->getAccelerations();
  velocities += 0.5 * dt * (accInit + accFinal);
  matter->setVelocities(velocities);
}

void Dynamics::run() {
  double sumT = 0.0, sumT2 = 0.0;

  setThermalVelocity();

  if (m_config.thermostat_kind != NONE) {
    QUILL_LOG_DEBUG(log,
                    "{} Running NVT molecular dynamics: {:8.2f} K for {} "
                    "steps ({:.4e} s)\n",
                    "[Dynamics]", temperature, m_config.steps,
                    1e-15 * m_config.time_step * m_config.timeUnit *
                        m_config.steps);
  } else {
    QUILL_LOG_DEBUG(log, "{} Running NVE molecular dynamics: {} steps\n",
                    "[Dynamics]", m_config.steps);
  }

  if (m_config.write_movies) {
    matter->matter2con("dynamics", false);
  }

  QUILL_LOG_DEBUG(log, "{} {:8} {:10} {:12} {:12} {:10}\n", "[Dynamics]",
                  "step", "KE", "PE", "TE", "kinT");

  for (long step = 0; step <= m_config.steps; step++) {
    oneStep();

    double kinE = matter->getKineticEnergy();
    double potE = matter->getPotentialEnergy();
    double kinT = 2.0 * kinE / nFreeCoords / kB;
    sumT += kinT;
    sumT2 += kinT * kinT;

    if (step % m_config.write_movies_interval == 0) {
      QUILL_LOG_DEBUG(log, "{} {} {} {} {} {}\n", "[Dynamics]", step, kinE,
                      potE, kinE + potE, kinT);
    }

    if (m_config.write_movies && (step % m_config.write_movies_interval == 0)) {
      matter->matter2con("dynamics", true);
    }
  }

  double avgT = sumT / static_cast<double>(m_config.steps);
  double varT = sumT2 / static_cast<double>(m_config.steps) - avgT * avgT;
  double stdT = std::sqrt(varT);
  QUILL_LOG_DEBUG(log,
                  "{} Temperature : Average = {:.2f} ; StdDev = {:.2f} ; "
                  "Factor = {:.2f}\n",
                  "[Dynamics]", avgT, stdT,
                  varT / avgT / avgT * nFreeCoords / 2.0);
}

void Dynamics::andersenCollision() {
  double alpha = m_config.andersen_alpha;
  double tCol = m_config.andersen_tcol;
  double pCol = 1.0 - std::exp(-m_config.time_step / tCol);

  AtomMatrix velocity = matter->getVelocities();
  auto mass = matter->getMasses();

  for (long i = 0; i < nAtoms; i++) {
    if (randomDouble() < pCol && !matter->getFixed(i)) {
      for (int j = 0; j < 3; j++) {
        double vOld = velocity(i, j);
        double vNew =
            std::sqrt(kB * temperature / mass[i]) * gaussRandom(0.0, 1.0);
        velocity(i, j) = std::sqrt(1.0 - alpha * alpha) * vOld + alpha * vNew;
      }
    }
  }
  matter->setVelocities(velocity);
}

void Dynamics::setThermalVelocity() {
  AtomMatrix velocity = matter->getVelocities();
  auto mass = matter->getMasses();

  for (long i = 0; i < nAtoms; i++) {
    if (!matter->getFixed(i)) {
      for (int j = 0; j < 3; j++) {
        velocity(i, j) =
            std::sqrt(kB * temperature / mass[i]) * gaussRandom(0.0, 1.0);
      }
    }
  }
  matter->setVelocities(velocity);
}

void Dynamics::rescaleVelocity() {
  AtomMatrix velocity = matter->getVelocities();
  double kinE = matter->getKineticEnergy();
  double kinT = 2.0 * kinE / nFreeCoords / kB;
  matter->setVelocities(velocity * std::sqrt(temperature / kinT));
}

/// Nose-Hoover chain thermostat (Martyna-Klein-Tuckerman algorithm).
/// Two chain variables (xi1, xi2) with velocities (vxi1, vxi2).
void Dynamics::noseHooverVerlet() {
  double dt2 = 0.5 * dt;
  double dt4 = 0.25 * dt;
  double dt8 = 0.125 * dt;
  double q1 = m_config.nose_mass;
  double q2 = q1;
  double Temp = kB * temperature;

  AtomMatrix vel = matter->getVelocities();
  AtomMatrix pos = matter->getPositions();
  double kinE = matter->getKineticEnergy();

  // Forward half-step for chain variables
  double g2 = (q1 * vxi1 * vxi1 - Temp);
  vxi2 += g2 * dt4;
  vxi1 *= std::exp(-vxi2 * dt8);
  double g1 = (2.0 * kinE - nFreeCoords * Temp) / q1;
  vxi1 += g1 * dt4;
  vxi1 *= std::exp(-vxi2 * dt8);
  xi1 += vxi1 * dt2;
  xi2 += vxi2 * dt2;
  double s = std::exp(-vxi1 * dt2);
  vel *= s;
  kinE *= s * s;
  vxi1 *= std::exp(-vxi2 * dt8);
  g1 = (2.0 * kinE - nFreeCoords * Temp) / q1;
  vxi1 += g1 * dt4;
  vxi1 *= std::exp(-vxi2 * dt8);
  g2 = (q1 * vxi1 * vxi1 - Temp) / q2;
  vxi2 += g2 * dt4;

  // Position + velocity Verlet step
  pos += vel * dt2;
  matter->setPositions(pos);
  AtomMatrix acc = matter->getAccelerations();
  vel += acc * dt;
  pos += vel * dt2;
  kinE = matter->getKineticEnergy();

  // Backward half-step for chain variables
  g2 = (q1 * vxi1 * vxi1 - Temp);
  vxi2 += g2 * dt4;
  vxi1 *= std::exp(-vxi2 * dt8);
  g1 = (2.0 * kinE - nFreeCoords * Temp) / q1;
  vxi1 += g1 * dt4;
  vxi1 *= std::exp(-vxi2 * dt8);
  xi1 += vxi1 * dt2;
  xi2 += vxi2 * dt2;
  s = std::exp(-vxi1 * dt2);
  vel *= s;
  kinE *= s * s;
  vxi1 *= std::exp(-vxi2 * dt8);
  g1 = (2.0 * kinE - nFreeCoords * Temp) / q1;
  vxi1 += g1 * dt4;
  vxi1 *= std::exp(-vxi2 * dt8);
  g2 = (q1 * vxi1 * vxi1 - Temp) / q2;
  vxi2 += g2 * dt4;

  matter->setPositions(pos);
  matter->setVelocities(vel);
}

/// Langevin dynamics (velocity-Verlet with friction and random forces).
void Dynamics::langevinVerlet() {
  double gamma = m_config.langevin_friction;
  AtomMatrix pos = matter->getPositions();
  AtomMatrix vel = matter->getVelocities();
  AtomMatrix acc = matter->getAccelerations();
  AtomMatrix noise = acc; // same shape
  auto mass = matter->getMasses();

  // Generate friction + stochastic forces
  AtomMatrix friction = -gamma * vel;
  for (long i = 0; i < nAtoms; i++) {
    if (!matter->getFixed(i)) {
      for (int j = 0; j < 3; j++) {
        noise(i, j) = std::sqrt(4.0 * gamma * kB * temperature / dt / mass[i]) *
                      gaussRandom(0.0, 1.0);
      }
    }
  }
  acc += friction + noise;

  vel += acc * 0.5 * dt;
  pos += vel * dt;
  matter->setPositions(pos);

  // Second half-step
  acc = matter->getAccelerations();
  friction = -gamma * vel;
  for (long i = 0; i < nAtoms; i++) {
    if (!matter->getFixed(i)) {
      for (int j = 0; j < 3; j++) {
        noise(i, j) = std::sqrt(4.0 * gamma * kB * temperature / dt / mass[i]) *
                      gaussRandom(0.0, 1.0);
      }
    }
  }
  acc += friction + noise;
  vel += 0.5 * dt * acc;
  matter->setVelocities(vel);
}
