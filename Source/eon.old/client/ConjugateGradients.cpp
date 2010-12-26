//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include "ConjugateGradients.h"
#include "HelperFunctions.h"

#include <cassert>
#include <cmath>

using namespace helper_functions;


ConjugateGradients::ConjugateGradients(Matter *matter, Parameters *parameters)
{
    initialize(matter, parameters);
    totalForceCalls = 0;
    outputLevel = 0;
}


ConjugateGradients::ConjugateGradients(Matter *matter,
                                       Parameters *parameters,
                                       Matrix<double, Eigen::Dynamic, 3> forces){
    initialize(matter, parameters);
    force = forces; 
}


void ConjugateGradients::initialize(Matter *matter_passed, Parameters *parameters_passed)
{
    // note that it is the pointer that is copied
    matter = matter_passed;
    parameters = parameters_passed;

    nAtoms =  matter->numberOfAtoms();

    direction.resize(nAtoms,3);
    directionOld.resize(nAtoms,3);
    directionNorm.resize(nAtoms,3);
    force.resize(nAtoms,3);
    forceOld.resize(nAtoms,3);

    direction.setZero();
    directionOld.setZero();
    directionNorm.setZero();

    force.setZero();
    forceOld.setZero();

    return;
}


ConjugateGradients::~ConjugateGradients(){

    // matter should not be deleted
    // parameters should not be deleted
    // Are pointers to objects outside the scope
 
    return;
}


void ConjugateGradients::oneStep(){
    long forceCallsTemp;
    double step;
    Matrix<double, Eigen::Dynamic, 3> pos;
    Matrix<double, Eigen::Dynamic, 3> posStep;
    Matrix<double, Eigen::Dynamic, 3> forceAfterStep;
    //----- Initialize end -----
    //std::cout<<"oneStep\n";

    forceCallsTemp = matter->getForceCalls();
    force = matter->getForces();
    pos = matter->getPositions();
    determineSearchDirection();

    // move system an infinitesimal step to determine the optimal step size along the search line
    posStep = pos + directionNorm * parameters->optFiniteDiffStep;
    matter->setPositions(posStep);
    forceAfterStep = matter->getForces();
 
    // move system optimal step
    step = stepSize(force, forceAfterStep, parameters->optMaxMove);
    //cout<<"Step: "<<step<<endl;
    pos += step*directionNorm;
    matter->setPositions(pos);

    forceCallsTemp = matter->getForceCalls()-forceCallsTemp;
    totalForceCalls += forceCallsTemp;

    return;
}


void ConjugateGradients::fullRelax(){
    bool converged = false;
    //----- Initialize end -----
    //std::cout<<"fullRelax\n";
    #ifndef NDEBUG
    static int run=0;
    ostringstream min;
    min << "min_" << run;
    matter->matter2xyz(min.str(), false);
    ++run;
    #endif
    int i=0;
    while(!converged and i < parameters->optMaxIterations) 
    {
        oneStep();
        converged = isItConverged(parameters->optConvergedForce);
        ++i;
        if (outputLevel > 0) {
            printf("step = %3d, max force = %8.5lf, energy: %10.4f\n", i, matter->maxForce(),
                   matter->getPotentialEnergy());
        }
        #ifndef NDEBUG
        if (i % 10 == 0)
            matter->matter2xyz(min.str(), true);
        #endif
    }
    return;
}


bool ConjugateGradients::isItConverged(double convergeCriterion){
    double diff=0;
 
    for(int i=0;i<nAtoms;i++)
    {
        diff = force.row(i).norm();
        if(convergeCriterion < diff)
        {
            break;
        }
    }
//    diff = length(force_,nFreeCoord_);
//    fprintf(stderr, "ConjugateGradients.isItConverged force magnitude: %f\n", diff);    
//std::cout<<diff<<"\n";
    return(diff < convergeCriterion);
}


void ConjugateGradients::setOutput(int level)
{
    outputLevel = level;
}


void ConjugateGradients::determineSearchDirection(){
    double a=0, b=0, gamma=0;
    //----- Initialize end -----
    //std::cout<<"determineSearchDirection\n";

    a = fabs((force.cwise() * forceOld).sum());
    b = forceOld.squaredNorm();
    if(a<0.5*b)
    {
        // Polak-Ribiere way to determine how much to mix in of old direction
        gamma = (force.cwise() * (force - forceOld)).sum()/b;
    }
    else
    {
        gamma = 0;
    }

    //std::cout << "forces:" << std::endl << force << std::endl;

    //cout<<"gamma: "<<gamma<<endl;
    direction = force + gamma*directionOld;
    //direction = direction.cwise() * matter->getFree();
    assert(direction.norm() != 0.0);
    directionNorm = direction;
    directionNorm.normalize();
    //cout<<"DirectionNorm: "<<directionNorm<<endl;
    directionOld = direction;
    forceOld = force;
    return;
}


double ConjugateGradients::stepSize(Matrix<double, Eigen::Dynamic, 3> forceBeforeStep, 
                                    Matrix<double, Eigen::Dynamic, 3> forceAfterStep,
                                    double maxStep){
    double projectedForce1;
    double projectedForce2;
    double step, curvature;
    //----- Initialize end -----
    //std::cout<<"stepSize\n";
    
    // Determine curvature
    projectedForce1 = (forceBeforeStep.cwise() * directionNorm).sum();
    projectedForce2 = (forceAfterStep.cwise() * directionNorm).sum();
    curvature = (projectedForce1-projectedForce2)/parameters->optFiniteDiffStep;
    
    if(curvature < 0)
        step = maxStep;
    else{
        step = projectedForce1/curvature;
        if(maxStep < fabs(step)){
            // Calculated is too large
            step = sign(step)*maxStep;
        }
    }
    return step;
}

// Specific functions when forces are modified 

Matrix<double, Eigen::Dynamic, 3> ConjugateGradients::makeInfinitesimalStepModifiedForces(Matrix<double, Eigen::Dynamic, 3> pos){
 
    determineSearchDirection();
    // Move system an infinitesimal step 
    // to determine the optimal step size along the search line
    return directionNorm * parameters->optFiniteDiffStep + pos;
}


Matrix<double, Eigen::Dynamic, 3> ConjugateGradients::getNewPosModifiedForces(
        Matrix<double, Eigen::Dynamic, 3> pos,
        Matrix<double, Eigen::Dynamic, 3> forceBeforeStep, 
        Matrix<double, Eigen::Dynamic, 3> forceAfterStep,
        double maxStep)
{
    double step;

    step = stepSize(forceBeforeStep, forceAfterStep, maxStep);

    // Move system
    return pos + directionNorm * step;
}


void ConjugateGradients::setForces(Matrix<double, Eigen::Dynamic, 3> forces){
    force = forces;
    return;
}
