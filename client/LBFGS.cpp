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

LBFGS::LBFGS(ObjectiveFunction *objfPassed, Parameters *parametersPassed)
{
    objf = objfPassed;
    parameters = parametersPassed;

    iteration = 0;

    //Shouldn't have a memory longer than the number of degrees of freedom.
    memory = min(objf->degreesOfFreedom(), (int)parameters->optLBFGSMemory);
}

LBFGS::~LBFGS()
{
    return;
}

VectorXd LBFGS::getStep(VectorXd f)
{
    double H0 = parameters->optLBFGSInverseCurvature;

    int loopmax = s.size();
    double a[loopmax];

    VectorXd q = -f;

    for (int i=loopmax-1;i>=0;i--) {
        a[i] = rho[i] * s[i].dot(q);
        q -= a[i] * y[i];
    }

    if (loopmax > 0 && parameters->optLBFGSAutoScale) {
        H0 = s[loopmax-1].dot(s[loopmax-1])/s[loopmax-1].dot(y[loopmax-1]);
        log_file("[LBFGS] H0: %.4e\n", H0); 
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

bool LBFGS::step(double maxMove)
{
    VectorXd r = objf->getPositions();
    VectorXd f = -objf->getGradient();

    if (iteration > 0) {
        update(r, rPrev, f, fPrev);
    }

    VectorXd d = getStep(f);
    double vd = d.normalized().dot(f.normalized());
    if (vd>1.0) vd=1.0;
    if (vd<-1.0) vd=-1.0;
    double angle = acos(vd) * (180.0 / M_PI);

    VectorXd dr;
    
    dr = helper_functions::maxAtomMotionAppliedV(d, maxMove);

    if (angle > 90.0) {
        log("LBFGS reset angle = %.4f\n", angle);
        s.erase(s.begin(), s.end());
        y.erase(y.begin(), y.end());
        rho.erase(rho.begin(), rho.end());
        d = getStep(f);
    }else if (dr != d) {
        //log("LBFGS reset, step too big\n");
        //s.erase(s.begin(), s.end());
        //y.erase(y.begin(), y.end());
        //rho.erase(rho.begin(), rho.end());
    }

    objf->setPositions(r+dr);

    rPrev = r;
    fPrev = f;

    iteration++;

    return objf->isConverged();
}


bool LBFGS::run(int maxSteps, double maxMove)
{
    while(!objf->isConverged() && iteration < maxSteps) {
        step(maxMove);
    }
    return objf->isConverged();
}
