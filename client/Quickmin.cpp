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
#include <cmath>

Quickmin::Quickmin(Matter *matter_passed, Parameters *parameters_passed)
{
    matter = matter_passed;
    parameters = parameters_passed;
    dt = parameters_passed->optTimeStep;
    velocity.resize(matter->numberOfAtoms(), 3);
    velocity.setZero();
    outputLevel = 0;
    totalForceCalls = 0;
}


Quickmin::~Quickmin()
{
    return;
}

void Quickmin::oneStep()
{
    AtomMatrix forces = matter->getForces();
    totalForceCalls++;
    if((velocity.cwise() * forces).sum() < 0)
    {
        dt *= 0.5;
        velocity.setZero();
    }
    else
    {
        dt *= 1.1;
    }
    velocity += forces * dt;
    AtomMatrix positions;
    positions = matter->getPositions();
    AtomMatrix step = velocity * dt;
    double maxmoved = 0.0;
    for(int i=0; i < matter->numberOfAtoms(); i++)
    {
        maxmoved = max(maxmoved, step.row(i).norm());
    }
    double scale = 1.0;
    if(maxmoved > parameters->optMaxMove)
    {
        scale = parameters->optMaxMove/maxmoved;
    }
    positions += step*scale;
    matter->setPositions(positions);  
}

void Quickmin::setOutput(int level)
{
    outputLevel = level;
}


long Quickmin::fullRelax()
{
    bool converged = false;
    ostringstream min;
    min << "min";
    if(parameters->writeMovies)
    {
        matter->matter2con(min.str(), false);
    }
    int i = 0;
    while(!converged)
    {
        if (i >= parameters->optMaxIterations) {
            return Minimizer::STATUS_MAX_ITERATIONS;
        }

        if (parameters->checkpointForceCalls != -1 &&
            Potential::fcalls >= parameters->checkpointForceCalls) {
            return Minimizer::STATUS_CHECKPOINT;
        }
        oneStep();
        converged = isItConverged(parameters->optConvergedForce);
        i++;
        if (outputLevel > 0) {
            printf("step = %3d, max force = %8.5lf, energy: %10.4f, dt: %8.5f\n", i,
                   matter->maxForce(), matter->getPotentialEnergy(), dt);
        }
        if(parameters->writeMovies)
        {
            matter->matter2con(min.str(), true);
        }
    }
    return Minimizer::STATUS_GOOD;
}


bool Quickmin::isItConverged(double convergeCriterion)
{
    AtomMatrix forces = matter->getForces();
    double diff = 0;
    for(int i=0;i<matter->numberOfAtoms();i++)
    {
        diff = forces.row(i).norm();
        if(convergeCriterion < diff)
        {
            break;
        }
    }
    return(diff < convergeCriterion);
}

