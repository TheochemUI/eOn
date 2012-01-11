//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include "FIRE.h"
#include "HelperFunctions.h"
#include <cmath>

FIRE::FIRE(ObjectiveFunction *objfPassed, Parameters *parametersPassed)
{
    objf = objfPassed;
    parameters = parametersPassed;
    dt = parametersPassed->optTimeStep;
    dt_max = parametersPassed->optMaxTimeStep;
    N_min = 5;
    N = 0;
    f_inc = 1.1;
    f_dec = 0.5;
    alpha_start = 0.1;
    alpha = alpha_start;
    f_a = 0.99;
    v0.resize(objf->degreesOfFreedom());
    v0.setZero();
    iteration = 0;
}

FIRE::~FIRE()
{
    return;
}

bool FIRE::step(double maxMove)
{
    // Check convergence.
    if(objf->isConverged()) {
        return true;
    }

    // Velocity Verlet
    VectorXd f0 = -objf->getGradient();
    VectorXd x0 = objf->getPositions();
    VectorXd dx = (v0 * dt) + (0.5 * f0 * dt * dt);
    dx = helper_functions::maxAtomMotionAppliedV(dx, parameters->optMaxMove);
    objf->setPositions(x0 + dx);
    VectorXd f1 = -objf->getGradient();
    VectorXd f1_unit = f1 / f1.norm();
    VectorXd v1 = v0 + ((f0 + f1) * 0.5) * dt;

    // FIRE
    double P = f1.dot(v1);
    v0 = (1 - alpha) * v1 + alpha * f1_unit * v1.norm();
    if(P >= 0) {
        N++;
        if (N > N_min) {
            dt = min(dt * f_inc, dt_max);
            alpha = alpha * f_a;
        }
    }
    else {
        dt = dt * f_dec;
        v0 = v0 * 0.0;
        alpha = alpha_start;
        N = 0;
    }

    iteration++;
    return objf->isConverged();
}

bool FIRE::run(int maxSteps, double maxMove)
{
    while(!objf->isConverged() && iteration < maxSteps) {
        step(maxMove);
    }
    return objf->isConverged();
}
