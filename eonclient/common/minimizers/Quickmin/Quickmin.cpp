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

Quickmin::Quickmin(Matter *matter, Parameters *parameters){
    // Note that it is the pointer that is copied.
    matter_ = matter;    
    parameters_ = parameters;
    
    nFreeCoord_ = 3*matter->numberOfFreeAtoms();
//    forces_ = new double[nFreeCoord_];
    tempListDouble_ = new double[nFreeCoord_];

};

Quickmin::~Quickmin(){

    // matter_ should not be deleted
    // parameters_ should not be deleted
    // forces_ should not be deleted
    // Are pointers to objects outside the scope
    
    delete [] tempListDouble_;
    return;
};


void Quickmin::oneStep(){
//    long forceCallsTemp;
    double *freeForces;
    freeForces = new double[nFreeCoord_];

    //----- Initialize end -----
    //std::cout<<"oneStep\n";

//    forceCallsTemp = matter_->getForceCalls();  
    matter_->getFreeForces(freeForces);
    oneStepPart1(freeForces);
    matter_->getFreeForces(freeForces);
    oneStepPart2(freeForces);

//    forceCallsTemp = matter_->getForceCalls()-forceCallsTemp;
//    parameters_->addForceCalls(forceCallsTemp);

    delete [] freeForces;
    return;
};

void Quickmin::oneStepPart1(double *freeForces){
    double *positions;
    double *velocity;
    positions = new double[nFreeCoord_];

    velocity = new double[nFreeCoord_];

    //----- Initialize end -----
    //std::cout<<"oneStepPart1\n";
    
    matter_->getFreeVelocities(velocity);  
    
    multiplyScalar(tempListDouble_, freeForces, 0.5*getTimeStep(), nFreeCoord_);
    add(velocity, tempListDouble_, velocity, nFreeCoord_);
    matter_->setFreeVelocities(velocity);
    
    matter_->getFreePositions(positions);       
    multiplyScalar(tempListDouble_, velocity, getTimeStep(), nFreeCoord_);
    add(positions, tempListDouble_, positions, nFreeCoord_);
    matter_->setFreePositions(positions);  

    delete [] positions;
    delete [] velocity;
    return;
};
    

void Quickmin::oneStepPart2(double *freeForces){
    double dotVelocityForces;
    double dotForcesForces;
    double *velocity;
    velocity = new double[nFreeCoord_];

    // Keep a copy of the pointer to the force array, is being used when it is 
    // decided if the calculation is converged
    forces_ = freeForces;
    //----- Initialize end -----
    //std::cout<<"oneStepPart2\n";
    
    matter_->getFreeVelocities(velocity); 
    multiplyScalar(tempListDouble_, freeForces, 0.5*getTimeStep(), nFreeCoord_);
    add(velocity, tempListDouble_, velocity, nFreeCoord_);
    
    dotVelocityForces = dot(velocity, freeForces, nFreeCoord_);
    // Zeroing all velocities if they are not orthogonal to the forces
    if(dotVelocityForces < 0)
        multiplyScalar(velocity, velocity, 0., nFreeCoord_);
    else{
        dotForcesForces = dot(freeForces, freeForces, nFreeCoord_);
        multiplyScalar(velocity, freeForces, 
                       dotVelocityForces/dotForcesForces, nFreeCoord_);
    }
    matter_->setFreeVelocities(velocity);      

    delete [] velocity;
    return;
};


void Quickmin::fullRelax(){
    bool converged = false;
    long forceCallsTemp;
    //----- Initialize end -----
    //std::cout<<"fullRelax\n";

    forceCallsTemp = matter_->getForceCalls();  
    while(!converged){
        oneStep();
        converged = isItConverged(parameters_->getConverged_Relax());
        std::cout<<matter_->potentialEnergy()<<"\n";
    }

    forceCallsTemp = matter_->getForceCalls()-forceCallsTemp;
    parameters_->addForceCalls(forceCallsTemp);
    return;
};


bool Quickmin::isItConverged(double convergeCriterion){
    double diff;
    
    for(int i=0;i<nFreeCoord_;i++){
        diff = fabs(forces_[i]);
        if(convergeCriterion < diff)
            break;
    }
    return(diff < convergeCriterion);
};


