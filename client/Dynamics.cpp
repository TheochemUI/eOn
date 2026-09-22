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
#include "eon/Dynamics.h"
#include "eon/EonLogger.h"

#include <cmath>
#include <stdexcept>

namespace eonc {

const char Dynamics::ANDERSEN[] = "andersen";
const char Dynamics::NOSE_HOOVER[] = "nose_hoover";
const char Dynamics::LANGEVIN[] = "langevin";
const char Dynamics::NONE[] = "none";

Dynamics::Dynamics(Matter *matter_in, const DynamicsConfig &config)
    : matter{matter_in},
      m_config{config} {
  if (!matter) {
    throw std::invalid_argument("Dynamics: null Matter");
  }
  dt = m_config.time_step;
  nAtoms = matter->numberOfAtoms();
  // Unfixed axes only. A partly fixed atom is not three degrees of freedom.
  nFreeCoords = 0;
  for (long i = 0; i < nAtoms; ++i) {
    for (int axis = 0; axis < 3; ++axis) {
      if (!matter->getFixed(i, axis)) {
        ++nFreeCoords;
      }
    }
  }
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
    const double kinT =
        (nFreeCoords > 0 && kB > 0.0) ? 2.0 * kinE / nFreeCoords / kB : 0.0;

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
    if (!eonc::io::io_ok(matter->matter2con("dynamics", false))) {
      QUILL_LOG_WARNING(log, "Failed to write dynamics movie header frame");
    }
  }

  QUILL_LOG_DEBUG(log, "{} {:8} {:10} {:12} {:12} {:10}\n", "[Dynamics]",
                  "step", "KE", "PE", "TE", "kinT");

  for (long step = 0; step < m_config.steps; step++) {
    oneStep();

    double kinE = matter->getKineticEnergy();
    double potE = matter->getPotentialEnergy();
    const double kinT =
        (nFreeCoords > 0 && kB > 0.0) ? 2.0 * kinE / nFreeCoords / kB : 0.0;
    sumT += kinT;
    sumT2 += kinT * kinT;

    if (step % m_config.write_movies_interval == 0) {
      QUILL_LOG_DEBUG(log, "{} {} {} {} {} {}\n", "[Dynamics]", step, kinE,
                      potE, kinE + potE, kinT);
    }

    if (m_config.write_movies && (step % m_config.write_movies_interval == 0)) {
      if (!eonc::io::io_ok(matter->matter2con("dynamics", true))) {
        QUILL_LOG_WARNING(log, "Failed to append dynamics movie frame");
      }
    }
  }

  const double nstat = static_cast<double>(std::max(m_config.steps, 1L));
  double avgT = sumT / nstat;
  double varT = sumT2 / nstat - avgT * avgT;
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
    if (eonc::rng::randomDouble() < pCol && !matter->getFixed(i)) {
      for (int j = 0; j < 3; j++) {
        double vOld = velocity(i, j);
        const double vNew = (mass[i] > 0.0 && kB > 0.0 && temperature > 0.0)
                                ? std::sqrt(kB * temperature / mass[i]) *
                                      eonc::rng::gaussRandom(0.0, 1.0)
                                : 0.0;
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
        velocity(i, j) = (mass[i] > 0.0 && kB > 0.0 && temperature > 0.0)
                             ? std::sqrt(kB * temperature / mass[i]) *
                                   eonc::rng::gaussRandom(0.0, 1.0)
                             : 0.0;
      }
    }
  }
  matter->setVelocities(velocity);
}

void Dynamics::rescaleVelocity() {
  AtomMatrix velocity = matter->getVelocities();
  double kinE = matter->getKineticEnergy();
  const double kinT =
      (nFreeCoords > 0 && kB > 0.0) ? 2.0 * kinE / nFreeCoords / kB : 0.0;
  if (!(kinT > 0.0) || !(temperature > 0.0)) {
    return;
  }
  matter->setVelocities(velocity * std::sqrt(temperature / kinT));
}

void Dynamics::nhcChainHalfStep(AtomMatrix &vel, double &kinE) {
  const double dt2 = 0.5 * dt;
  const double dt4 = 0.25 * dt;
  const double dt8 = 0.125 * dt;
  const double q1 = m_config.nose_mass;
  const double q2 = q1;
  const double Temp = kB * temperature;
  if (!(q1 > 0.0)) {
    throw std::invalid_argument("thermostat.nose_mass must be positive");
  }

  // Martyna, Klein, Tuckerman JCP 97, 2635 (1992): G2 = (Q1 v_ξ1² − kT) / Q2.
  auto g2 = [&]() { return (q1 * vxi1 * vxi1 - Temp) / q2; };
  auto g1 = [&]() { return (2.0 * kinE - nFreeCoords * Temp) / q1; };

  vxi2 += g2() * dt4;
  vxi1 *= std::exp(-vxi2 * dt8);
  vxi1 += g1() * dt4;
  vxi1 *= std::exp(-vxi2 * dt8);
  xi1 += vxi1 * dt2;
  xi2 += vxi2 * dt2;
  const double s = std::exp(-vxi1 * dt2);
  vel *= s;
  kinE *= s * s;
  vxi1 *= std::exp(-vxi2 * dt8);
  vxi1 += g1() * dt4;
  vxi1 *= std::exp(-vxi2 * dt8);
  vxi2 += g2() * dt4;
}

/// Nose-Hoover chain thermostat (Martyna-Klein-Tuckerman algorithm).
/// Two chain variables (xi1, xi2) with velocities (vxi1, vxi2).
void Dynamics::noseHooverVerlet() {
  const double dt2 = 0.5 * dt;
  AtomMatrix vel = matter->getVelocities();
  AtomMatrix pos = matter->getPositions();
  double kinE = matter->getKineticEnergy();

  nhcChainHalfStep(vel, kinE);

  pos += vel * dt2;
  matter->setPositions(pos);
  AtomMatrix acc = matter->getAccelerations();
  vel += acc * dt;
  pos += vel * dt2;
  matter->setPositions(pos);
  matter->setVelocities(vel);
  kinE = matter->getKineticEnergy();

  nhcChainHalfStep(vel, kinE);

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
                      eonc::rng::gaussRandom(0.0, 1.0);
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
                      eonc::rng::gaussRandom(0.0, 1.0);
      }
    }
  }
  acc += friction + noise;
  vel += 0.5 * dt * acc;
  matter->setVelocities(vel);
}

} // namespace eonc
