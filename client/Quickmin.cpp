
/*
 *===============================================
 *  Created by Andreas Pedersen on 4/17/07.
 *-----------------------------------------------
 *  Modified. Name, Date and a small description!
 *
 *-----------------------------------------------
 *  Todo:
 *
 *===============================================
 */

#include "Quickmin.h"

using namespace helper_functions;
using namespace constants;


Quickmin::Quickmin(Matter *matter, Parameters *parameters)
{
    matter_ = matter;    
    parameters_ = parameters;
    dtScale_ = 1.0;
    
    nFreeCoord_ = 3*matter->numberOfFreeAtoms();
    tempListDouble_ = new double[nFreeCoord_];

};


Quickmin::~Quickmin()
{
    /* matter_, parameters_, and forces_ should not be deleted. They are pointers to objects outside the scope.*/
    delete [] tempListDouble_;
    return;
};


void Quickmin::oneStep()
{
    double *freeForces;
    freeForces = new double[nFreeCoord_];
    matter_->getFreeForces(freeForces);
    oneStepPart1(freeForces);
    matter_->getFreeForces(freeForces);
    oneStepPart2(freeForces);
    delete [] freeForces;
    return;
};


void Quickmin::oneStepPart1(double *freeForces)
{
    double *positions;
    double *velocity;
    positions = new double[nFreeCoord_];
    velocity = new double[nFreeCoord_];
    matter_->getFreeVelocities(velocity);  
    multiplyScalar(tempListDouble_, freeForces, 0.5 * getTimeStep() * dtScale_, nFreeCoord_);
    add(velocity, tempListDouble_, velocity, nFreeCoord_);
    matter_->setFreeVelocities(velocity);
    matter_->getFreePositions(positions);       
    multiplyScalar(tempListDouble_, velocity, getTimeStep() * dtScale_, nFreeCoord_);
    add(positions, tempListDouble_, positions, nFreeCoord_);
    matter_->setFreePositions(positions);  
    delete [] positions;
    delete [] velocity;
    return;
};
    

void Quickmin::oneStepPart2(double *freeForces)
{
    double dotVelocityForces;
    double dotForcesForces;
    double *velocity;
    velocity = new double[nFreeCoord_];
    forces_ = freeForces;
    matter_->getFreeVelocities(velocity); 
    multiplyScalar(tempListDouble_, freeForces, 0.5*getTimeStep() * dtScale_, nFreeCoord_);
    add(velocity, tempListDouble_, velocity, nFreeCoord_);
    dotVelocityForces = dot(velocity, freeForces, nFreeCoord_);
    // Zeroing all velocities if they are not orthogonal to the forces
    if(dotVelocityForces < 0)
    {
        multiplyScalar(velocity, velocity, 0., nFreeCoord_);
        dtScale_ *= 0.99;
    }    
    else
    {
        dotForcesForces = dot(freeForces, freeForces, nFreeCoord_);
        multiplyScalar(velocity, freeForces, dotVelocityForces/dotForcesForces, nFreeCoord_);
//        dtScale_ *= 1.01;
    }
    matter_->setFreeVelocities(velocity);      
    delete [] velocity;
    return;
};


void Quickmin::fullRelax()
{
    bool converged = false;
    long forceCallsTemp;
    forceCallsTemp = matter_->getForceCalls();  
    while(!converged)
    {
        oneStep();
        converged = isItConverged(parameters_->getConverged_Relax());
    }
    forceCallsTemp = matter_->getForceCalls()-forceCallsTemp;
    parameters_->addForceCalls(forceCallsTemp);
    return;
};


bool Quickmin::isItConverged(double convergeCriterion)
{
    double diff;
    for(int i=0;i<nFreeCoord_;i++)
    {
        diff = fabs(forces_[i]);
        if(convergeCriterion < diff)
        {
            break;
        }
    }
    return(diff < convergeCriterion);
};


