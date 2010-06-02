/*
 *===============================================
 *  Created by Andreas Pedersen on 10/11/06.
 *-----------------------------------------------
 *  Modified. Name, Date and a small description!
 *
 *-----------------------------------------------
 *  Todo:
 *
 *-----------------------------------------------
 *  Heavily inspired of codes by:
 *      Graeme Henkelmann
 *      Andri Arnaldsson
 *      Roar Olsen
 *===============================================
 */
#include <cassert>
#include "ConjugateGradients.h"

using namespace helper_functions;
using namespace constants;


ConjugateGradients::ConjugateGradients(Matter *matter, Parameters *parameters){
    initialize(matter, parameters);
};


// Constructor to be used when modified forces are used
ConjugateGradients::ConjugateGradients(Matter *matter, 
                                       Parameters *parameters, 
                                       double *forces){
    initialize(matter, parameters);
    
    for(int i=0; i<nFreeCoord_; i++)
        force_[i] = forces[i];    
};


void ConjugateGradients::initialize(Matter *matter, Parameters *parameters)
{
    // Note that it is the pointer that is copied.
    matter_ = matter;    
    parameters_ = parameters;
    
    nFreeCoord_ = 3 * matter->numberOfFreeAtoms();
    
    direction_ = new double[nFreeCoord_];
    directionOld_ = new double[nFreeCoord_];
    directionNorm_ = new double[nFreeCoord_];

    force_ = new double[nFreeCoord_];
    forceOld_ = new double[nFreeCoord_];
    
    // There should be space for both free and frozen atoms
    tempListDouble_ = new double[3*matter->numberOfAtoms()];

    for(int i=0;i<nFreeCoord_;i++)
        direction_[i]=directionOld_[i]=directionNorm_[i]=force_[i]=forceOld_[i]=0;

    return;
};


ConjugateGradients::~ConjugateGradients(){

    // matter_ should not be deleted
    // parameters_ should not be deleted
    // Are pointers to objects outside the scope
    
    delete [] tempListDouble_;
    delete [] direction_;
    delete [] directionOld_;
    delete [] directionNorm_;
    delete [] force_;
    delete [] forceOld_;
    return;
};


void ConjugateGradients::oneStep(){
    long forceCallsTemp;
    double step;
    double *pos;
    double *posStep;
    double *forceAfterStep;
    pos = new double[nFreeCoord_];
    posStep = new double[nFreeCoord_];
    forceAfterStep = new double[nFreeCoord_];
    //----- Initialize end -----
    //std::cout<<"oneStep\n";

    forceCallsTemp = matter_->getForceCalls();    
    matter_->getFreeForces(force_);
    assert(length(force_, nFreeCoord_) != 0.0);
    matter_->getFreePositions(pos);  
    determineSearchDirection();
    // Move system an infinitisimal step 
    // to determine the optimal step size along the search line
    multiplyScalar(tempListDouble_, directionNorm_, 
                   getCurvatureStep(), nFreeCoord_);
    
    add(posStep, tempListDouble_, pos, nFreeCoord_);
    matter_->setFreePositions(posStep);
    matter_->getFreeForces(forceAfterStep);
    
    // Move system optimal step
    step = stepSize(force_, forceAfterStep, getMaxMoveFullRelax());
    multiplyScalar(tempListDouble_, directionNorm_, step, nFreeCoord_);
    add(pos, tempListDouble_, pos, nFreeCoord_);    
    matter_->setFreePositions(pos);

    forceCallsTemp = matter_->getForceCalls()-forceCallsTemp;
    parameters_->addForceCalls(forceCallsTemp);
    
    delete [] pos;
    delete [] posStep;
    delete [] forceAfterStep;
    return;
};


void ConjugateGradients::fullRelax(){
    bool converged = false;
    //----- Initialize end -----
    //std::cout<<"fullRelax\n";
    int i=0;
    while(!converged and i < parameters_->getMaximumIterations()) {
        oneStep();
        converged = isItConverged(parameters_->getConverged_Relax());
        ++i;
        #ifndef NDEBUG
        printf("min = %d\n", i);
        #endif
    }
    #ifndef NDEBUG
    matter_->matter2con("minimum");
    #endif
    return;
};


bool ConjugateGradients::isItConverged(double convergeCriterion){
    double diff;
    
//    for(int i=0;i<nFreeCoord_;i++){
//        diff = fabs(force_[i]);//-forceOld_[i]);
//        if(convergeCriterion < diff)
//            break;
//    }
    diff = length(force_,nFreeCoord_);
//    fprintf(stderr, "ConjugateGradients.isItConverged force magnitude: %f\n", diff);    
//std::cout<<diff<<"\n";

    return(diff < convergeCriterion);
};


void ConjugateGradients::determineSearchDirection(){
    assert(length(force_, nFreeCoord_) != 0.0);
    double a=0, b=0, gamma=0;
    //----- Initialize end -----
    //std::cout<<"determineSearchDirection\n";

    a = fabs(dot(force_,forceOld_,nFreeCoord_));
    b = dot(forceOld_,forceOld_,nFreeCoord_);
    if(a<0.5*b){
        subtract(tempListDouble_, force_, forceOld_, nFreeCoord_);
        //Polak-Ribiere way to determine how much to mix in of old direction
        gamma = dot(force_, tempListDouble_, nFreeCoord_)/b;  
    }
    else
        gamma = 0;
    
    multiplyScalar(tempListDouble_, directionOld_, gamma, nFreeCoord_);
    add(direction_, force_, tempListDouble_, nFreeCoord_);
    assert(length(direction_, nFreeCoord_) != 0.0);
    copyRightIntoLeft(directionNorm_, direction_, nFreeCoord_);
    normalize(directionNorm_, nFreeCoord_);
    
    copyRightIntoLeft(directionOld_, direction_, nFreeCoord_);
    copyRightIntoLeft(forceOld_, force_, nFreeCoord_);
    return;
};


double ConjugateGradients::stepSize(double *forceBeforeStep, 
                                    double *forceAfterStep,
                                    double maxStep){
    double projectedForce1;
    double projectedForce2;
    double step, curvature;
    //----- Initialize end -----
    //std::cout<<"stepSize\n";
    
    // Determine curvature
    projectedForce1 = dot(forceBeforeStep,directionNorm_,nFreeCoord_);
    projectedForce2 = dot(forceAfterStep,directionNorm_,nFreeCoord_);
    curvature = (projectedForce1-projectedForce2)/getCurvatureStep();
    
    if(curvature < 0)
        step = maxStep;
    else{
        step = projectedForce1/curvature;
        if(maxStep < fabs(step)){
            // Too big step was calculated
            step = sign(step)*maxStep;
        }
    }
    return step;
};

// Specific functions when forces are modified 

void ConjugateGradients::makeInfinitisimalStepModifiedForces(double *posStep, 
                                                             double *pos){
    
    determineSearchDirection();
    
    // Move system an infinitisimal step 
    // to determine the optimal step size along the search line
    multiplyScalar(tempListDouble_,
                   directionNorm_,
                   getCurvatureStep(),
                   nFreeCoord_);
    add(posStep, tempListDouble_, pos, nFreeCoord_);
    return;
};


void ConjugateGradients::getNewPosModifiedForces(double *pos,
                                                 double *forceBeforeStep, 
                                                 double *forceAfterStep,
                                                 double maxStep){
    double step;

    step = stepSize(forceBeforeStep, forceAfterStep, maxStep);
    
    // Move system
    multiplyScalar(tempListDouble_, directionNorm_, step, nFreeCoord_);
    add(pos, tempListDouble_, pos, nFreeCoord_);
    return;
};

    
void ConjugateGradients::setFreeAtomForcesModifiedForces(double *forces){
    for(int i=0; i<nFreeCoord_; i++)
        force_[i] = forces[i];
    return;
};
