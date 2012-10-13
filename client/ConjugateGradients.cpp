//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include "ConjugateGradients.h"
#include "Log.h"
#include <cassert>
#include <cmath>

ConjugateGradients::ConjugateGradients(ObjectiveFunction *objf_in, Parameters *parameters_in)
{
    objf = objf_in;
    parameters = parameters_in;
    forceOld = objf->getPositions() * 0.0;
    directionOld = objf->getPositions() * 0.0;
}

VectorXd ConjugateGradients::getStep()
{
    double a=0, b=0, gamma=0;
    a = fabs(force.dot(forceOld));
    b = forceOld.squaredNorm();
    if (a < 0.5 * b) {
        // Polak-Ribiere way to determine how much to mix in of old direction
        gamma = force.dot(force - forceOld) / b;
    } else {
        gamma = 0;
    }
    direction = force + gamma * directionOld;
    directionNorm = direction;
    directionNorm.normalize();
    directionOld = direction;
    forceOld = force;
    return direction;
}

bool ConjugateGradients::step(double maxMove)
{
    VectorXd pos;
    VectorXd posStep;
    VectorXd forceAfterStep;

    force = -objf->getGradient();
    pos = objf->getPositions();
    getStep();

    // move system an infinitesimal step to determine the optimal step size along the search line
    posStep = pos + directionNorm * parameters->finiteDifference;
    objf->setPositions(posStep);
    forceAfterStep = -objf->getGradient(true);

    // Determine curvature
    double projectedForce1 = force.dot(directionNorm);
    double projectedForce2 = forceAfterStep.dot(directionNorm);
    double curvature = (projectedForce1 - projectedForce2) / parameters->finiteDifference;

    double stepSize = 100.0;

    if(curvature > 0.0)
    {
        stepSize = projectedForce1 / curvature;
    }
    if (!parameters->optCGNoOvershooting)
    {
        pos += helper_functions::maxAtomMotionAppliedV(stepSize * directionNorm, maxMove);
        objf->setPositions(pos);
    }
    else
    {
        // negative if product of the projected forces before and after the step are in opposite directions
        double passedMinimum = -1.;
        double forceChange = 0.;
        while (passedMinimum < 0.) {
            posStep = pos + helper_functions::maxAtomMotionAppliedV(stepSize * directionNorm, maxMove);
            objf->setPositions(posStep);
            forceAfterStep = -objf->getGradient(true);
            projectedForce2 = forceAfterStep.dot(directionNorm);
            
            passedMinimum = projectedForce1 * projectedForce2;
            if (passedMinimum < 0.)
            {
                forceChange = (projectedForce1 - projectedForce2);
                stepSize = (projectedForce1 / forceChange) * stepSize;
                
                // knockout old search direction
                directionOld = objf->getPositions() * 0.0;
            }
        }
    }
    return objf->isConverged();
}


bool ConjugateGradients::run(int maxIterations, double maxMove)
{
    int iterations = 0;
    while(!objf->isConverged() && iterations <= maxIterations)
    {
        step(maxMove);
        iterations++;
    }
    return objf->isConverged();
}


ConjugateGradients::~ConjugateGradients()
{
    return;
}

