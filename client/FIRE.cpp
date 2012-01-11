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

VectorXd FIRE::getStep()
{
    VectorXd a0 = -objf->getGradient();
    VectorXd x0 = objf->getPositions();
    VectorXd x1 = x0 + (v0 * dt) + (0.5 * a0 * dt * dt);
    objf->setPositions(x1);
    VectorXd a1 = -objf->getGradient();
    objf->setPositions(x0);
    VectorXd a1_unit = a1 / a1.norm();
    VectorXd v1 = v0 + ((a0 + a1) * 0.5) * dt;
    double P = a1.dot(v1);
    v0 = (1 - alpha) * v1 + alpha * a1_unit * v1.norm();
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
    return x1 - x0;
}

bool FIRE::step(double maxMove)
{
    if(objf->isConverged()) {
        return true;
    }
    VectorXd d = getStep();
    VectorXd dr = helper_functions::maxAtomMotionAppliedV(d, parameters->optMaxMove);
    objf->setPositions(objf->getPositions() + dr);  
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
