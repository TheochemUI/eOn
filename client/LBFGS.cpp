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

    if (iteration > 0 && parameters->optLBFGSAutoScale) {
        VectorXd dr = objf->getPositions() - rPrev;
        VectorXd dg = -f+fPrev;
        H0 = dr.dot(dr)/dr.dot(dg);
    }
    if ((iteration == 0 || H0 < 0) && parameters->optLBFGSAutoScale) {
        //calculate H0 via finite difference
        VectorXd r = objf->getPositions();
        objf->setPositions(r+parameters->finiteDifference*f.normalized());
        VectorXd dg = objf->getGradient(true)+f;
        H0 = dg.dot(f.normalized())/parameters->finiteDifference;
        H0 = 1.0/H0;
        objf->setPositions(r);
        if (H0 > 0) {
            log_file("[LBFGS] H0 calculated via FD: %.4e\n", H0); 
        }else{
            log_file("[LBFGS] H0 calculated via FD: %.4e, sd instead\n", H0); 
            iteration = 0;
            return 100*f.normalized();
        }

    }else if (parameters->optLBFGSAutoScale) {
        log_file("[LBFGS] H0: %.4e\n", H0); 
    }

    int loopmax = s.size();
    double a[loopmax];

    VectorXd q = -f;

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

    //y0 is the change in the gradient, not the force
    VectorXd y0 = f0 - f1;

    s.push_back(s0);
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
        double C = (r-rPrev).dot(fPrev-f)/(r-rPrev).dot(r-rPrev);
        if (C<0) {
            log("[LBFGS] Curvature: %.4f\n",C);
            s.clear();
            y.clear();
            rho.clear();
            iteration = 0;
        }
    }


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

    bool reset;
    reset = false;
    if (angle > 90.0 && parameters->optLBFGSAngleReset) {
        log("LBFGS reset, angle, %.4f\n", angle);
        reset = true;
    }

    if (dr != d && parameters->optLBFGSDistanceReset) {
        log("LBFGS reset, step too big, %.4f\n", d.norm());
        reset = true;
    }

    if (reset) {
        s.clear();
        y.clear();
        rho.clear();
        d = getStep(f);
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
