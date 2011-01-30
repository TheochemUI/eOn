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
#include "Constants.h"
#include "Parameters.h"

#include <cmath>
using namespace helper_functions;


Quickmin::Quickmin(Matter *matter_passed, Parameters *parameters_passed)
{
    matter = matter_passed;    
    parameters = parameters_passed;
    dt = parameters_passed->optTimeStep;
    velocity.resize(matter->numberOfAtoms(), 3);
    velocity.setZero();
    outputLevel = 0;    
}


Quickmin::~Quickmin()
{
    /* matter_, parameters_, and forces_ should not be deleted. They are pointers to objects outside the scope.*/
    return;
}

void Quickmin::oneStep()
{
    Matrix<double, Eigen::Dynamic, 3> forces = matter->getForces();
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
    Matrix<double, Eigen::Dynamic, 3> positions;
    positions = matter->getPositions();
    Matrix<double, Eigen::Dynamic, 3> step = velocity * dt;
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


void Quickmin::fullRelax()
{
    bool converged = false;
    long forceCallsTemp;
    forceCallsTemp = matter->getForceCalls();  
    int i = 0;
    while(!converged)
    {
        oneStep();
        converged = isItConverged(parameters->optConvergedForce);
        i++;
        if (outputLevel > 0) {
            printf("step = %3d, max force = %8.5lf, energy: %10.4f, dt: %8.5f\n", i, matter->maxForce(),
                   matter->getPotentialEnergy(), dt);
        }
    }
    forceCallsTemp = matter->getForceCalls()-forceCallsTemp;
    return;
}


bool Quickmin::isItConverged(double convergeCriterion)
{
    Matrix<double, Eigen::Dynamic, 3> forces = matter->getForces();
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

