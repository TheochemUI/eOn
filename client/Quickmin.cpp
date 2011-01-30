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
    dtScale = 1.0;
    
    nAtoms = matter->numberOfAtoms();
    velocity.resize(nAtoms, 3);
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
    Matrix<double, Eigen::Dynamic, 3> forces;
    forces = matter->getForces();
    oneStepPart1(forces);
    forces = matter->getForces();
    oneStepPart2(forces);
    return;
}

void Quickmin::setOutput(int level)
{
    outputLevel = level;
}


void Quickmin::oneStepPart1(Matrix<double, Eigen::Dynamic, 3> force)
{
    Matrix<double, Eigen::Dynamic, 3> positions;

    velocity += force * .5 * parameters->optTimeStep * dtScale;
    velocity = velocity.cwise() * matter->getFree();

    positions = matter->getPositions();
    Matrix<double, Eigen::Dynamic, 3> update = velocity * parameters->optTimeStep * dtScale;
    double maxmoved = 0.0;
    for(int i=0; i < matter->numberOfAtoms(); i++)
    {
        maxmoved = max(maxmoved, update.row(i).norm());
    }
    double scale = 1.0;
    if(maxmoved > parameters->optMaxMove)
    {
        scale = parameters->optMaxMove/maxmoved;
    }
    positions += update*scale;
    matter->setPositions(positions);  
}


void Quickmin::oneStepPart2(Matrix<double, Eigen::Dynamic, 3> force)
{
    double dotVelocityForces;
    double dotForcesForces;
    forces = force;
    velocity += force * 0.5 * parameters->optTimeStep * dtScale;
    velocity = velocity.cwise() * matter->getFree();
    dotVelocityForces = (velocity.cwise() * force).sum();
    // Zeroing all velocities if they are not orthogonal to the forces
    if(dotVelocityForces < 0)
    {
        velocity.setZero();
        dtScale *= 0.99;
    }    
    else
    {
        dotForcesForces = force.squaredNorm();
        velocity = force * dotVelocityForces/dotForcesForces;
//        dtScale_ *= 1.01;
    }
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
            printf("step = %3d, max force = %8.5lf, energy: %10.4f\n", i, matter->maxForce(),
                   matter->getPotentialEnergy());
        }
    }
    forceCallsTemp = matter->getForceCalls()-forceCallsTemp;
    return;
}


bool Quickmin::isItConverged(double convergeCriterion)
{
    double diff = 0;
    for(int i=0;i<nAtoms;i++)
    {
        diff = forces.row(i).norm();
        if(convergeCriterion < diff)
        {
            break;
        }
    }
    return(diff < convergeCriterion);
}

