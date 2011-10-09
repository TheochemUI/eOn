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
#include "Log.h"
#include <cassert>
#include <cmath>

ConjugateGradients::ConjugateGradients(Matter *matter, Parameters *parameters)
{
    initialize(matter, parameters);
    totalForceCalls = 0;
    outputLevel = 0;
}


ConjugateGradients::ConjugateGradients(Matter *matter, Parameters *parameters, AtomMatrix forces)
{
    initialize(matter, parameters);
    totalForceCalls = 0;
    outputLevel = 0;
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


ConjugateGradients::~ConjugateGradients()
{
    // matter should not be deleted
    // parameters should not be deleted
    // Are pointers to objects outside the scope

    return;
}


void ConjugateGradients::oneStep()
{
    long forceCallsTemp;
    AtomMatrix pos;
    AtomMatrix posStep;
    AtomMatrix forceAfterStep;
    //----- Initialize end -----
    //std::cout<<"oneStep\n";

    forceCallsTemp = matter->getForceCalls();
    force = matter->getForces();
    pos = matter->getPositions();
    determineSearchDirection();
    // move system an infinitesimal step to determine the optimal step size along the search line
    posStep = pos + directionNorm * parameters->finiteDifference;
    matter->setPositions(posStep);
    forceAfterStep = matter->getForces();

    //cout<<"Step: "<<step<<endl;
    pos += getStep(force, forceAfterStep, parameters->optMaxMove);
    matter->setPositions(pos);

    forceCallsTemp = matter->getForceCalls()-forceCallsTemp;
    totalForceCalls += forceCallsTemp;

    return;
}


long ConjugateGradients::fullRelax()
{
    //----- Initialize end -----
    //std::cout<<"fullRelax\n";
    ostringstream min;
    min << "min";
    if(parameters->writeMovies)
    {
        if (parameters->checkpoint) {
            matter->matter2con(min.str(), true);
        }else{
            matter->matter2con(min.str(), false);
        }
    }
    int i = 0;
    force = matter->getForces();
    bool converged = isItConverged(parameters->optConvergedForce);
    while(!converged)
    {
        if (i >= parameters->optMaxIterations) {
            return Minimizer::STATUS_MAX_ITERATIONS;
        }

        oneStep();
        converged = isItConverged(parameters->optConvergedForce);
        ++i;

        if (!parameters->quiet) {
            log("step = %3d, max force = %10.7lf, energy: %10.7f\n", 
                   i, matter->maxForce(), matter->getPotentialEnergy());
        }
        if (parameters->writeMovies) {
            matter->matter2con(min.str(), true);
        }
        if (parameters->checkpoint) {
            matter->matter2con("reactant_checkpoint.con");
        }

    }
    return Minimizer::STATUS_GOOD;
}


bool ConjugateGradients::isItConverged(double convergeCriterion)
{
    double diff = 0;

    for(int i=0; i<nAtoms; i++)
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
    if(a < 0.5*b)
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


AtomMatrix ConjugateGradients::getStep(AtomMatrix forceBeforeStep, AtomMatrix forceAfterStep, double maxStep, bool saddleSearch)
{
    double projectedForce1;
    double projectedForce2;
    double stepSize, curvature;
    //----- Initialize end -----
    //std::cout<<"stepSize\n";

    // Determine curvature
    projectedForce1 = (forceBeforeStep.cwise() * directionNorm).sum();
    projectedForce2 = (forceAfterStep.cwise() * directionNorm).sum();
    curvature = (projectedForce1-projectedForce2)/parameters->finiteDifference;

    // new
    Matrix<double, Eigen::Dynamic, 3> directionTemp;
    directionTemp = directionNorm;
    // new
//    assert(fabs(saddleConfinePositive) < 1.0);
    if(curvature < 0.0)
    {
        stepSize = 100.0;
        if(parameters->saddleConfinePositive) {
            // new
            // Zero out the directions with the least displacement
            int moved = 0;
            double minMoved = parameters->saddleConfinePositiveMinMove;
            while(moved == 0){
                moved = 0;
                directionTemp = directionNorm;
                for (int i=0; i<3*nAtoms; i++) {
                    if (fabs(directionTemp[i]) < minMoved)
                        directionTemp[i]=0;
                    else{
                        moved = moved + 1;
                    }
                }
                minMoved *= parameters->saddleConfinePositiveScaleRatio;
            }

            if (moved > parameters->saddleConfinePositiveMaxActiveAtoms){
                printf("crash, too many atoms are moving (increase confine_positive_min_move? )\n");
                exit(1);
            }
            std::cout<<moved<<" components are moving"<<std::endl;
        }

    }
    else
    {
        stepSize = projectedForce1/curvature;
    }

//    directionTemp.normalize();
//    return stepSize * directionTemp;
//    return helper_functions::maxAtomMotionApplied(stepSize * directionTemp, maxStep);
      return helper_functions::maxAtomMotionApplied(stepSize * directionTemp, maxStep);
}

/*
    if(curvature < 0.0)
        stepSize = maxStep;
    else{
        stepSize = projectedForce1/curvature;
        if(maxStep < fabs(stepSize)){
            // Calculated is too large
            #ifndef NDEBUG
                cout<<"CG exceeded max step"<<endl;
            #endif
            double sign = (stepSize > 0.0 ) ? 1.0 : -1.0;
            stepSize = sign*maxStep;
        }
    }
    return stepSize * directionNorm;
}
*/
// Specific functions when forces are modified 

AtomMatrix ConjugateGradients::makeInfinitesimalStepModifiedForces(AtomMatrix pos){
 
    determineSearchDirection();
    // Move system an infinitesimal step 
    // to determine the optimal step size along the search line
    return directionNorm * parameters->finiteDifference + pos;
}


AtomMatrix ConjugateGradients::getNewPosModifiedForces(
        AtomMatrix pos,
        AtomMatrix forceBeforeStep,
        AtomMatrix forceAfterStep,
        double maxStep,
        bool saddleSearch)
{
    // Move system
    return pos + getStep(forceBeforeStep, forceAfterStep, maxStep, saddleSearch);
}


void ConjugateGradients::setForces(AtomMatrix forces){
    force = forces;
    return;
}
