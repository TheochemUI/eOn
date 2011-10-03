//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include "QuickminBox.h"
#include <cmath>

QuickminBox::QuickminBox(Matter *matter_passed, Parameters *parameters_passed)
{
    matter = matter_passed;
    parameters = parameters_passed;
    dt = parameters_passed->optTimeStep;
    velocity.resize(matter->numberOfAtoms(), 3);
    velocity.setZero();
    boxXv = boxYv = boxZv = boxXf = boxYf = boxZf = 0.0;
    outputLevel = 0;
}

QuickminBox::~QuickminBox()
{
    /* matter_, parameters_, and forces_ should not be deleted. 
    They are pointers to objects outside the scope.*/
    return;
}

void QuickminBox::oneStep()
{
    double dR = parameters->optFiniteDist;
    Vector3d boxv = matter->getBoundary(0);
    matter->setBoundary(0, boxv.normalized() * (boxv.norm() + dR));
    double eX = matter->getPotentialEnergy();
    matter->setBoundary(0, boxv);
    boxv = matter->getBoundary(1);
    matter->setBoundary(1, boxv.normalized() * (boxv.norm() + dR));
    double eY = matter->getPotentialEnergy();
    matter->setBoundary(1, boxv);
    boxv = matter->getBoundary(2);
    matter->setBoundary(2, boxv.normalized() * (boxv.norm() + dR));
    double eZ = matter->getPotentialEnergy();
    matter->setBoundary(2, boxv);
    double e0 = matter->getPotentialEnergy();
    boxXf = -(eX - e0)/dR;
    boxYf = -(eY - e0)/dR;
    boxZf = -(eZ - e0)/dR;

    AtomMatrix forces = matter->getForces();
    if((velocity.cwise() * forces).sum() < 0 || (boxXv * boxXf + boxYv * boxYf + boxZv * boxZf) < 0.0)
    {
        dt *= 0.5;
        velocity.setZero();
        boxXv = boxYv = boxZv = 0.0;
    }
    else
    {
        dt *= 1.1;
    }

    boxXv += boxXf * dt;
    boxYv += boxYf * dt;
    boxZv += boxZf * dt;
    double stepX = min(parameters->optMaxMove, max(-parameters->optMaxMove, boxXv * dt));
    double stepY = min(parameters->optMaxMove, max(-parameters->optMaxMove, boxYv * dt));
    double stepZ = min(parameters->optMaxMove, max(-parameters->optMaxMove, boxZv * dt));
    boxv = matter->getBoundary(0);
    matter->setBoundary(0, boxv.normalized() * (boxv.norm() + stepX));
    boxv = matter->getBoundary(1);
    matter->setBoundary(1, boxv.normalized() * (boxv.norm() + stepY));
    boxv = matter->getBoundary(2);
    matter->setBoundary(2, boxv.normalized() * (boxv.norm() + stepZ));

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
    positions += step * scale;
    matter->setPositions(positions);  
}

void QuickminBox::setOutput(int level)
{
    outputLevel = level;
}


long QuickminBox::fullRelax()
{
    bool converged = false;
    long forceCallsTemp;
    forceCallsTemp = matter->getForceCalls();  
    int i = 0;
    while(!converged)
    {
        if (i >= parameters->optMaxIterations) {
            return Minimizer::STATUS_MAX_ITERATIONS;
        }

        oneStep();
        converged = isItConverged(parameters->optConvergedForce);
        i++;
        if (outputLevel > 0) {
            printf("step = %3d, max force = %8.5lf, energy: %10.4f, dt: %8.5f, boxForce: % 8.5f, box: % 3.3f % 3.3f % 3.3f\n", i, 
                   matter->maxForce(), matter->getPotentialEnergy(), dt, sqrt(boxXf*boxXf + boxYf*boxYf + boxZf*boxZf), 
                   matter->getBoundary(0).norm(), matter->getBoundary(1).norm(), matter->getBoundary(2).norm());
        }
    }
    forceCallsTemp = matter->getForceCalls()-forceCallsTemp;
    return Minimizer::STATUS_GOOD;
}


bool QuickminBox::isItConverged(double convergeCriterion)
{
    AtomMatrix forces = matter->getForces();
    double diff = 0;
    if(sqrt(boxXf*boxXf + boxYf*boxYf + boxZf*boxZf) > convergeCriterion)
    {
        return false;
    }
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

