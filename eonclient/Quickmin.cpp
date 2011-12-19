//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include "Quickmin.h"
#include "HelperFunctions.h"
#include <cmath>

Quickmin::Quickmin(ObjectiveFunction *objfPassed, Parameters *parametersPassed)
{
    objf = objfPassed;
    parameters = parametersPassed;
    dt = parametersPassed->optTimeStep;
    velocity.resize(objf->degreesOfFreedom());
    velocity.setZero();
    iteration = 0;
}

Quickmin::~Quickmin()
{
    return;
}

VectorXd Quickmin::getStep()
{
    VectorXd forces = -objf->getGradient();

    if (velocity.dot(forces) < 0) {
        if (parameters->optVariableTimeStep) dt *= 0.5;
        velocity.setZero();
    } else {
        if (parameters->optVariableTimeStep) dt *= 1.1;
    }
    velocity += forces * dt;
    return velocity * dt;
}

bool Quickmin::step(double maxMove)
{
    VectorXd d = getStep();
    VectorXd dr = helper_functions::maxAtomMotionAppliedV(d, parameters->optMaxMove);

    VectorXd positions = objf->getPositions();
    positions += dr;
    objf->setPositions(positions);  
    iteration++;
    return objf->isConverged();
}

bool Quickmin::run(int maxSteps, double maxMove)
{
    while(!objf->isConverged() && iteration < maxSteps) {
        step(maxMove);
    }
    return objf->isConverged();
}
