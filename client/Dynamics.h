#ifndef DYNAMICS_H
#define DYNAMICS_H

#include "Optimizer.h"
#include "Matter.h"
#include "HelperFunctions.h"
#include "Parameters.h"

#include "Eigen.h"

class Dynamics {

public:

    Dynamics(Matter *matter, Parameters *parameters);

    ~Dynamics();

    void setTemperature(double temperature);
    void oneStep(int stepNumber=-1);
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
};

#endif
