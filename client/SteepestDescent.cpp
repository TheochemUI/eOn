//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

//Based on the SteepestDescent minimizer written in ASE.

#include "SteepestDescent.h"
#include "Log.h"
#include <cassert>
#include <cmath>

SteepestDescent::SteepestDescent(ObjectiveFunction *objfPassed, Parameters *parametersPassed)
{
    objf = objfPassed;
    parameters = parametersPassed;

    iteration = 0;
}

SteepestDescent::~SteepestDescent()
{
    return;
}

bool SteepestDescent::step(double maxMove)
{
    VectorXd r = objf->getPositions();
    VectorXd f = -objf->getGradient();

    VectorXd dr;
    double alpha = parameters->optSDAlpha;
    if (parameters->optSDTwoPoint == true && iteration > 0) {
        VectorXd dx = r-rPrev;
        VectorXd dg = -f+fPrev;
        alpha = dx.dot(dx)/dx.dot(dg);
        if (alpha < 0) {
            alpha = parameters->optSDAlpha;
        }
        log_file("[SD] alpha: %.4e\n", alpha);
    }

    dr = alpha*f;
    dr = helper_functions::maxAtomMotionAppliedV(dr, maxMove);

    objf->setPositions(r+dr);

    rPrev = r;
    fPrev = f;

    iteration++;

    return objf->isConverged();
}


bool SteepestDescent::run(int maxSteps, double maxMove)
{
    while(!objf->isConverged() && iteration < maxSteps) {
        step(maxMove);
    }
    return objf->isConverged();
}
