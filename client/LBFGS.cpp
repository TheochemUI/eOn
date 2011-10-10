//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

//Based on the LBFGS minimizer written in ASE.

#include "LBFGS.h"
#include "Log.h"
#include <cassert>
#include <cmath>

LBFGS::LBFGS(Matter *matter, Parameters *parametersPassed)
{
    totalForceCalls = 0;
    outputLevel = 0;
    parameters = parametersPassed;
    objf = new MatterObjectiveFunction(matter, parameters);
    ePrev = 0;
    iteration = 0;

    //Shouldn't have a memory longer than the number of degrees of freedom.
    memory = min(objf->degreesOfFreedom(), (int)parameters->optLBFGSMemory);
}


LBFGS::LBFGS(Matter *matter, Parameters *parameters, AtomMatrix forces)
{
}


LBFGS::~LBFGS()
{
    delete objf;
    return;
}

VectorXd LBFGS::getDescentDirection()
{
    double H0 = 1./10.;

    int loopmax = s.size();
    double a[loopmax];

    VectorXd q = objf->getGradient();

    for (int i=loopmax-1;i>=0;i--) {
        a[i] = rho[i] * s[i].dot(q);
        q -= a[i] * y[i];
    }

    VectorXd z = H0 * q;

    for (int i=0;i<loopmax;i++) {
        double b = rho[i] * y[i].dot(z);
        z += s[i] * (a[i] - b);
    }

    VectorXd d = -z;

    return d;
}


void LBFGS::update(VectorXd r1, VectorXd r0, VectorXd f1, VectorXd f0)
{
    VectorXd s0 = r1 - r0;
    s.push_back(s0);

    //y0 is the change in the gradient, not the force
    VectorXd y0 = f0 - f1;
    y.push_back(y0);

    rho.push_back(1.0/(s0.dot(y0)));

    if ((int)s.size() > memory) {
        s.erase(s.begin());                                                                          
        y.erase(y.begin());
        rho.erase(rho.begin());
    }
}

void LBFGS::oneStep()
{
    VectorXd r = objf->getPositions();
    VectorXd f = -objf->getGradient();

    if (iteration > 0) {
        update(r, rPrev, f, fPrev);
    }

    VectorXd d = getDescentDirection();
    double vd = d.normalized().dot(f.normalized());
    if (vd>1.0) vd=1.0;
    double angle = acos(vd) * (180.0 / M_PI);
    //printf("\nangle: %.3f\n", angle);

    if (angle > 70.0) {
        int keep = min((int)s.size(), 1);
        s.erase(s.begin(), s.end()-keep);
        y.erase(y.begin(), y.end()-keep);
        rho.erase(rho.begin(), rho.end()-keep);
        d = getDescentDirection();
    }

    VectorXd dr;
    dr = helper_functions::maxAtomMotionAppliedV(d, parameters->optMaxMove);

    objf->setPositions(r+dr);

    rPrev = r;
    fPrev = f;

    iteration++;

    return;
}


long LBFGS::fullRelax()
{
    while(!objf->converged())
    {
        if (iteration >= parameters->optMaxIterations) {
            return Minimizer::STATUS_MAX_ITERATIONS;
        }

        oneStep();

        if (!parameters->quiet) {
            double e = objf->getEnergy();
            log("step = %3d, max_force: %10.4f energy: %10.4f de: %10.2e memory: %i\n", 
                iteration, objf->convergence(), e, e-ePrev, s.size());
            ePrev = objf->getEnergy();
        }

    }
    return Minimizer::STATUS_GOOD;
}


bool LBFGS::isItConverged(double convergeCriterion)
{
    VectorXd gradient = objf->getGradient();
    for (int i=0; i<objf->degreesOfFreedom();i++) {
        if (convergeCriterion < gradient(i)) {
            return false;
        }
    }
    return true;
}


void LBFGS::setOutput(int level)
{
    outputLevel = level;
}
