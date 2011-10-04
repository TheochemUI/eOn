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
    memory = parameters->optLBFGSMemory;

    return;
}


LBFGS::~LBFGS()
{
    return;
}


void LBFGS::oneStep()
{
    VectorXd r = matter->getPositionsFreeV();
    VectorXd f = matter->getForcesFreeV();

    double H0 = 1./10.;

    if (iteration > 0) {
        VectorXd s0 = r - r0;
        s.push_back(s0);
        VectorXd y0 = f0 - f;
        y.push_back(y0);
        rho.push_back(1.0/(s0.dot(y0)));
    }

    if (iteration > memory) {
        s.pop_back();
        y.pop_back();
        rho.pop_back();
    }

    int loopmax = min(iteration, memory);
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

    VectorXd dr = -z;


    dr = helper_functions::maxAtomMotionAppliedV(dr, parameters->optMaxMove);

    matter->setPositionsFreeV(r+dr);

    r0 = r;
    f0 = f;

    force = matter->getForces();

    iteration++;

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
            log("step = %3d, max force = %10.7lf, energy: %10.7f\n", 
                   i, matter->maxForce(), matter->getPotentialEnergy());
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
