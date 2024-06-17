#include "Dynamics.h"
#include "Eigen.h"
#include <math.h>

using namespace helper_functions;

const char Dynamics::ANDERSEN[] = "andersen";
const char Dynamics::NOSE_HOOVER[] = "nose_hoover";
const char Dynamics::LANGEVIN[] = "langevin";
const char Dynamics::NONE[] = "none";

Dynamics::Dynamics(Matter *matter_in, Parameters *parameters_in) {
  matter = matter_in;
  parameters = parameters_in;
  dt = parameters->md.timeStep;
  nAtoms = matter->numberOfAtoms();
  nFreeCoords = matter->numberOfFreeAtoms() * 3;
  temperature = parameters->main.temperature;
  kB = parameters->kB;
  vxi1 = vxi2 = xi1 = xi2 = 0.0; // NoseHoover variables
  log = spdlog::get("combi");
}

Dynamics::~Dynamics() { return; }

void Dynamics::setTemperature(double temperature_in) {
  temperature = temperature_in;
}

void Dynamics::oneStep(int stepNumber) {
  // TODO(rg): These should be named enums
  if (parameters->thermostat.kind == ANDERSEN) {
    andersenCollision();
    velocityVerlet();
  } else if (parameters->thermostat.kind == NOSE_HOOVER) {
    noseHooverVerlet();
  } else if (parameters->thermostat.kind == LANGEVIN) {
    langevinVerlet();
  } else if (parameters->thermostat.kind == NONE) {
    velocityVerlet();
  }

  if (stepNumber != -1) {
    if (stepNumber == 1) {
      SPDLOG_LOGGER_DEBUG(log, "{} {:8s} {:10s} {:12s} {:12s} {:10s}\n",
                          "[Dynamics]", "Step", "KE", "PE", "TE", "KinT");
    }
    AtomMatrix velocity;
    double potE, kinE, kinT;
    velocity = matter->getVelocities();
    kinE = matter->getKineticEnergy();
    potE = matter->getPotentialEnergy();
    kinT = (2.0 * kinE / nFreeCoords / kB);

    if (stepNumber % parameters->debug.writeMoviesInterval == 0) {
      SPDLOG_LOGGER_DEBUG(
          log, "{} {:8} {:10.4} {:12.4} {:12.4} {:10.2}\n", "[Dynamics]",
          stepNumber, kinE, potE, kinE + potE, kinT);
    }
  }
}

void Dynamics::velocityVerlet() {
  AtomMatrix positions = matter->getPositions();
  AtomMatrix velocities = matter->getVelocities();
  AtomMatrix accelerationsInitial = matter->getAccelerations();

  positions += (dt * velocities) + (0.5 * dt * dt * accelerationsInitial);
  matter->setPositions(positions);

  AtomMatrix accelerationsFinal = matter->getAccelerations();

  velocities += 0.5 * dt * (accelerationsInitial + accelerationsFinal);
  matter->setVelocities(velocities);
}

void Dynamics::run() {
  AtomMatrix velocity;
  double potE, kinE, kinT;
  double sumT = 0, sumT2 = 0, avgT, varT, stdT;

  setThermalVelocity();

  if (parameters->thermostat.kind != NONE) {
    SPDLOG_LOGGER_DEBUG(log,
                        "{} Running NVT molecular dynamics: {:8.2lf} K for {} "
                        "steps ({:.4e} s)\n",
                        "[Dynamics]", temperature, parameters->md.steps,
                        1e-15 * parameters->md.timeStep * parameters->timeUnit *
                            parameters->md.steps);
  } else {
    SPDLOG_LOGGER_DEBUG(log, "{} Running NVE molecular dynamics: {} steps\n",
                        "[Dynamics]", parameters->md.steps);
  }

  if (parameters->debug.writeMovies == true) {
    matter->matter2con("dynamics", false);
  }

  SPDLOG_LOGGER_DEBUG(log, "%s %8s %10s %12s %12s %10s\n", "[Dynamics]"s,
                      "step", "KE", "PE", "TE", "kinT");

  for (long step = 0; step <= parameters->md.steps; step++) {
    oneStep();

    velocity = matter->getVelocities();
    kinE = matter->getKineticEnergy();
    potE = matter->getPotentialEnergy();
    kinT = (2.0 * kinE / nFreeCoords / kB);
    sumT += kinT;
    sumT2 += kinT * kinT;

    if (step % parameters->debug.writeMoviesInterval == 0) {
      SPDLOG_LOGGER_DEBUG(log, "{} {:8} {:10} {:12} {:12} {:10}\n",
                          "[Dynamics]", step, kinE, potE, kinE + potE, kinT);
    }

    if ((parameters->debug.writeMovies == true) &&
        (step % parameters->debug.writeMoviesInterval == 0)) {
      matter->matter2con("dynamics", true);
    }
  }
  avgT = sumT / double(parameters->md.steps);
  varT = sumT2 / double(parameters->md.steps) - avgT * avgT;
  stdT = sqrt(varT);
  SPDLOG_LOGGER_DEBUG(log,
                      "{} Temperature : Average = {:.2lf} ; StdDev = {:.2lf} ; "
                      "Factor = {:.2lf}\n",
                      "[Dynamics]", avgT, stdT,
                      varT / avgT / avgT * nFreeCoords / 2.0);
}

void Dynamics::andersenCollision() {
  double alpha, tCol, pCol;
  double vNew, vOld;
  VectorType mass;
  AtomMatrix velocity;

  alpha = parameters->thermostat.andersenAlpha; // collision strength
  tCol = parameters->thermostat.andersenTcol; // average time between collisions, in
                                         // unit of fs
  pCol = 1.0 - exp(-parameters->md.timeStep / tCol);

  velocity = matter->getVelocities();
  mass = matter->getMasses();

  for (long i = 0; i < nAtoms; i++) {
    if ((randomDouble() < pCol) && (!matter->getFixed(i))) {
      for (int j = 0; j < 3; j++) {
        vOld = velocity(i, j);
        vNew = sqrt(kB * temperature / mass[i]) * gaussRandom(0.0, 1.0);
        velocity(i, j) = sqrt(1.0 - alpha * alpha) * vOld + alpha * vNew;
      }
    }
  }
  matter->setVelocities(velocity);
}

void Dynamics::setThermalVelocity() {
  AtomMatrix velocity = matter->getVelocities();
  VectorType mass = matter->getMasses();

  for (long i = 0; i < nAtoms; i++) {
    if (!matter->getFixed(i)) {
      for (int j = 0; j < 3; j++) {
        velocity(i, j) =
            sqrt(kB * temperature / mass[i]) * gaussRandom(0.0, 1.0);
      }
    }
  }
  matter->setVelocities(velocity);
}

void Dynamics::rescaleVelocity() {
  AtomMatrix velocity = matter->getVelocities();
  double kinE = matter->getKineticEnergy();
  double kinT = (2.0 * kinE / nFreeCoords / kB);
  matter->setVelocities(velocity * sqrt(temperature / kinT));
}

void Dynamics::noseHooverVerlet() {
  AtomMatrix vel, pos, acc;
  double q1, q2, g1, g2, s, dt2, dt4, dt8;
  double kinE, Temp;

  dt2 = 0.5 * dt;
  dt4 = 0.25 * dt;
  dt8 = 0.125 * dt;
  q1 = q2 = parameters->thermostat.noseMass;
  g1 = g2 = 0.0;
  Temp = kB * temperature; // imposed temperature

  vel = matter->getVelocities();
  pos = matter->getPositions();
  acc = matter->getAccelerations();

  kinE = matter->getKineticEnergy();

  g2 = (q1 * vxi1 * vxi1 - Temp);
  vxi2 += g2 * dt4;
  vxi1 *= exp(-vxi2 * dt8);
  g1 = (2.0 * kinE - nFreeCoords * Temp) / q1;
  vxi1 += g1 * dt4;
  vxi1 *= exp(-vxi2 * dt8);
  xi1 += vxi1 * dt2;
  xi2 += vxi2 * dt2;
  s = exp(-vxi1 * dt2);
  vel *= s;
  kinE *= s * s;
  vxi1 *= exp(-vxi2 * dt8);
  g1 = (2.0 * kinE - nFreeCoords * Temp) / q1;
  vxi1 += g1 * dt4;
  vxi1 *= exp(-vxi2 * dt8);
  g2 = (q1 * vxi1 * vxi1 - Temp) / q2;
  vxi2 += g2 * dt4;

  pos += vel * dt2;
  matter->setPositions(pos);

  acc = matter->getAccelerations();
  vel += acc * dt;
  pos += vel * dt2;
  kinE = matter->getKineticEnergy();

  g2 = (q1 * vxi1 * vxi1 - Temp);
  vxi2 += g2 * dt4;
  vxi1 *= exp(-vxi2 * dt8);
  g1 = (2.0 * kinE - nFreeCoords * Temp) / q1;
  vxi1 += g1 * dt4;
  vxi1 *= exp(-vxi2 * dt8);
  xi1 += vxi1 * dt2;
  xi2 += vxi2 * dt2;
  s = exp(-vxi1 * dt2);
  vel *= s;
  kinE *= s * s;
  vxi1 *= exp(-vxi2 * dt8);
  g1 = (2.0 * kinE - nFreeCoords * Temp) / q1;
  vxi1 += g1 * dt4;
  vxi1 *= exp(-vxi2 * dt8);
  g2 = (q1 * vxi1 * vxi1 - Temp) / q2;
  vxi2 += g2 * dt4;

  matter->setPositions(pos);
  matter->setVelocities(vel);
}

void Dynamics::langevinVerlet() {
  VectorType mass;
  AtomMatrix pos;
  AtomMatrix vel;
  AtomMatrix acc;
  AtomMatrix friction;
  AtomMatrix noise;
  double gamma;

  gamma = parameters->thermostat.langevinFriction;
  pos = matter->getPositions();
  vel = matter->getVelocities();

  acc = matter->getAccelerations();
  noise = acc;
  mass = matter->getMasses();

  friction = -gamma * vel;
  for (long i = 0; i < nAtoms; i++) {
    if (!matter->getFixed(i)) {
      for (int j = 0; j < 3; j++) {
        noise(i, j) = sqrt(4.0 * gamma * kB * temperature / dt / mass[i]) *
                      gaussRandom(0.0, 1.0);
      }
    }
  }
  acc = acc + friction + noise;

  vel += acc * 0.5 * dt; // calculate velocites v(n+1/2)
  pos += vel * dt;
  matter->setPositions(pos); // update positions x(n+1)

  acc = matter->getAccelerations();
  friction = -gamma * vel;
  for (long i = 0; i < nAtoms; i++) {
    if (!matter->getFixed(i)) {
      for (int j = 0; j < 3; j++) {
        noise(i, j) = sqrt(4.0 * gamma * kB * temperature / dt / mass[i]) *
                      gaussRandom(0.0, 1.0);
      }
    }
  }
  acc = acc + friction + noise;
  vel += 0.5 * dt * acc;
  matter->setVelocities(vel); // calculate velocities v(n+1)
}
