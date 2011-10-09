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

LBFGS::LBFGS(Matter *matter, Parameters *parameters)
{
    initialize(matter, parameters);
    totalForceCalls = 0;
    outputLevel = 0;
}


LBFGS::LBFGS(Matter *matter, Parameters *parameters, AtomMatrix forces)
{
    initialize(matter, parameters);
    totalForceCalls = 0;
    outputLevel = 0;
}


void LBFGS::initialize(Matter *matter_passed, Parameters *parameters_passed)
{
    // note that it is the pointer that is copied
    matter = matter_passed;
    parameters = parameters_passed;

    iteration = 0;
    degreesOfFreedom = 3*matter->numberOfFreeAtoms();

    //Shouldn't have a memory longer than the number of degrees of freedom.
    memory = min(degreesOfFreedom, (int)parameters->optLBFGSMemory);

    return;
}


LBFGS::~LBFGS()
{
    return;
}

void LBFGS::update(VectorXd r1, VectorXd r0, VectorXd f1, VectorXd f0)
{
    if (iteration > 0) {
        VectorXd s0 = r1 - r0;
        s.push_back(s0);

        //y0 is the change in the gradient, not the force
        VectorXd y0 = f0 - f1;
        y.push_back(y0);

        rho.push_back(1.0/(y0.dot(s0)));
    }

    if (iteration > memory) {
        s.erase(s.begin());
        y.erase(y.begin());
        rho.erase(rho.begin());
    }
}

VectorXd LBFGS::getDescentDirection()
{
    double H0 = 1./10.;

    int loopmax = s.size();
    double a[loopmax];

    VectorXd q = -matter->getForcesFreeV();

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

void LBFGS::lineSearch(VectorXd d)
{
    double alphak = 1.0;
    double sigma = 1e-2;
    double scale = 0.90;

    double e0 = matter->getPotentialEnergy();
    VectorXd g0 = -matter->getForcesFreeV();
    VectorXd r0 = matter->getPositionsFreeV();

    matter->setPositionsFreeV(r0 + alphak * d);
    double e = matter->getPotentialEnergy();

    int i=0;
    while (e > e0+sigma*alphak*g0.dot(d) && i < 4) {
        alphak = alphak*scale;
        matter->setPositionsFreeV(r0 + alphak*d);
        e = matter->getPotentialEnergy();
        i++;
    }
    if (i>0) printf("alphak: %12.4f e: %12.4f e0: %12.4f\n", alphak, e, e0);
    return;
}


void LBFGS::oneStep()
{
    VectorXd r = matter->getPositionsFreeV();
    VectorXd f = matter->getForcesFreeV();

    update(r, rPrev, f, fPrev);

    VectorXd d = getDescentDirection();
    VectorXd dr = helper_functions::maxAtomMotionAppliedV(d, parameters->optMaxMove);
    matter->setPositionsFreeV(r+dr);

    rPrev = r;
    fPrev = f;

    iteration++;

    force = matter->getForces();

    return;
}


long LBFGS::fullRelax()
{
    string min = "min";
    if(parameters->writeMovies)
    {
        if (parameters->checkpoint) {
            matter->matter2con(min.c_str(), true);
        }else{
            matter->matter2con(min.c_str(), false);
        }
    }

    force = matter->getForces();
    bool converged = isItConverged(parameters->optConvergedForce);

    int i=0;
    while(!converged)
    {
        if (i >= parameters->optMaxIterations) {
            return Minimizer::STATUS_MAX_ITERATIONS;
        }

        oneStep();

        converged = isItConverged(parameters->optConvergedForce);
        i++;

        if (!parameters->quiet) {
            double e = matter->getPotentialEnergy();
            if (iteration > 0) {
                log("step = %3d, max force = %10.7lf, energy: %10.7f de: %10.2e\n", 
                    i, matter->maxForce(), e, e-ePrev);
            }else{
                log("step = %3d, max force = %10.7lf, energy: %10.7f\n", 
                    i, matter->maxForce(), e);
            }
            ePrev = matter->getPotentialEnergy();
        }
        if (parameters->writeMovies) {
            matter->matter2con(min.c_str(), true);
        }
        if (parameters->checkpoint) {
            matter->matter2con("reactant_checkpoint.con");
        }

    }
    return Minimizer::STATUS_GOOD;
}


bool LBFGS::isItConverged(double convergeCriterion)
{
    double diff = 0;

    for(int i=0; i<matter->numberOfAtoms(); i++)
    {
        diff = force.row(i).norm();
        if(convergeCriterion < diff)
        {
            break;
        }
    }

    return(diff < convergeCriterion);
}


void LBFGS::setOutput(int level)
{
    outputLevel = level;
}


AtomMatrix LBFGS::getStep(AtomMatrix forceBeforeStep, AtomMatrix forceAfterStep, double maxStep, bool saddleSearch)
{
    AtomMatrix tmp;
    return tmp;
}

AtomMatrix LBFGS::makeInfinitesimalStepModifiedForces(AtomMatrix pos){
    AtomMatrix tmp;
    return tmp;
}


AtomMatrix LBFGS::getNewPosModifiedForces(
        AtomMatrix pos,
        AtomMatrix forceBeforeStep,
        AtomMatrix forceAfterStep,
        double maxStep,
        bool saddleSearch)
{
    return pos + getStep(forceBeforeStep, forceAfterStep, maxStep, saddleSearch);
}


void LBFGS::setForces(AtomMatrix forces){
    force = forces;
    return;
}
