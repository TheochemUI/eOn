#ifndef DYNAMICS_H
#define DYNAMICS_H

#include "HelperFunctions.h"
#include "Matter.h"
#include "Optimizer.h"
#include "Parameters.h"

#include "Eigen.h"

class Dynamics {

public:
  static const char ANDERSEN[];
  static const char NOSE_HOOVER[];
  static const char LANGEVIN[];
  static const char NONE[];

  Dynamics(Matter *matter, Parameters *parameters);

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
  long nAtoms, nFreeCoords;

  Matter *matter;
  Parameters *parameters;

  double dt;
  double kB;
  double temperature;
  double vxi1, vxi2, xi1, xi2;
  std::shared_ptr<spdlog::logger> log;
};

#endif
