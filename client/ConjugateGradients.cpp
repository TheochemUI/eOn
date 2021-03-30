#include "ConjugateGradients.h"
#include "Log.h"
#include <cassert>
#include <cmath>


ConjugateGradients::ConjugateGradients(ObjectiveFunction *objf_in, Parameters *parameters_in)
{
    objf = objf_in;
    parameters = parameters_in;
    forceOld = objf->getPositions() * 0.0;
    directionOld = objf->getPositions() * 0.0;
    cg_i = 0;
}

VectorXd ConjugateGradients::getStep()
{
    double a=0, b=0, gamma=0;
    a = fabs(force.dot(forceOld));
    b = forceOld.squaredNorm();
    if (a < 0.5 * b) {
        // Polak-Ribiere way to determine how much to mix in of old direction
        gamma = force.dot(force - forceOld) / b;
    } else {
        gamma = 0;
    }
    direction = force + gamma * directionOld;
    directionNorm = direction;
    directionNorm.normalize();
    directionOld = direction;
    forceOld = force;
    
    // Only if value for max nr of iteration before reset
    if ((parameters->optCGMaxIterBeforeReset > 0) and
        (parameters->optCGMaxIterBeforeReset <= cg_i ))
        //or gamma == 0))
    {
        cg_i = 0;
        forceOld = objf->getPositions() * 0.0;
        directionOld = objf->getPositions() * 0.0;
//        std::cout<<"reset\n";
    }
    cg_i += 1;
    
    return direction;
}


int ConjugateGradients::step(double maxMove)
{
    bool converged;
    if (parameters->optCGLineSearch)
    {
        converged = line_search(maxMove);
    }
    else
    {
        converged = single_step(maxMove);
    }
    if(converged) return 1;
    return 0;
}


int ConjugateGradients::line_search(double maxMove)
{
    VectorXd pos;
    VectorXd posStep;
    VectorXd forceBeforeStep;
    double stepSize;
    double projectedForce;
    double projectedForceBeforeStep;
    double curvature;
    
    forceBeforeStep = -objf->getGradient();
    force = forceBeforeStep;
    getStep();
    pos = objf->getPositions();
    projectedForceBeforeStep = force.dot(directionNorm);
    
    // move system an infinitesimal step to determine the optimal step size along the search line
    posStep = pos + directionNorm * parameters->finiteDifference;
    objf->setPositions(posStep);
    force = -objf->getGradient(true);
    projectedForce = force.dot(directionNorm);
    stepSize = parameters->finiteDifference;
    
    int line_i = 0;
    do
    {
        // Determine curvature from last step (Secant method)
        curvature = fabs((projectedForceBeforeStep - projectedForce) / stepSize);
        stepSize = projectedForce / curvature;
        //stepSize = projectedForceBeforeStep / curvature;
        
        if (maxMove < fabs(stepSize))
        {
            // first part get the sign of stepSize
            stepSize = ((stepSize > 0) - (stepSize < 0)) * maxMove;
        }
        
        forceBeforeStep = force;
        projectedForceBeforeStep = projectedForce;
        
        pos += stepSize * directionNorm;
        objf->setPositions(pos);
        force = -objf->getGradient();
        projectedForce = force.dot(directionNorm);
        
        line_i += 1;
        
    // Line search considered converged based in the ratio between the projected force and the norm of the true force 
    } while (parameters->optCGLineConverged < fabs(projectedForce) / (sqrt(force.dot(force)+parameters->optCGLineConverged))
             and (line_i < parameters->optCGLineSearchMaxIter));
//    return objf->isConverged();
    if(objf->isConverged()) return 1;
    return 0;
}


int ConjugateGradients::single_step(double maxMove)
{
    VectorXd pos;
    VectorXd posStep;
    VectorXd forceAfterStep;
    
    force = -objf->getGradient();
    pos = objf->getPositions();
    getStep();
    
    // move system an infinitesimal step to determine the optimal step size along the search line
    posStep = pos + directionNorm * parameters->finiteDifference;
    objf->setPositions(posStep);
    forceAfterStep = -objf->getGradient(true);
    
    // Determine curvature
    double projectedForce1 = force.dot(directionNorm);
    double projectedForce2 = forceAfterStep.dot(directionNorm);
    double curvature = (projectedForce1 - projectedForce2) / parameters->finiteDifference;
    
    double stepSize = maxMove;
    
    if(curvature > 0.0){
        stepSize = projectedForce1 / curvature;
    }

    if(parameters->saddleBowlBreakout and maxMove < 0.0){
        stepSize = -maxMove;
        maxMove = -maxMove;
    }
    
    if (!parameters->optCGNoOvershooting)
    {
        if(parameters->saddleBowlBreakout){
             // max displacement is based on system not single atom
             pos += helper_functions::maxMotionAppliedV(stepSize * directionNorm, maxMove);
        }
        else{
            pos += helper_functions::maxAtomMotionAppliedV(stepSize * directionNorm, maxMove);
        }
        objf->setPositions(pos);
    }
    else
    {
        // negative if product of the projected forces before and after the step are in opposite directions
        double passedMinimum = -1.;
        double forceChange = 0.;
        while (passedMinimum < 0. and (0.1 * fabs(projectedForce1) < fabs(projectedForce2)))
        {
            posStep = pos + helper_functions::maxAtomMotionAppliedV(stepSize * directionNorm, maxMove);
            objf->setPositions(posStep);
            forceAfterStep = -objf->getGradient(true);
            projectedForce2 = forceAfterStep.dot(directionNorm);
            
            passedMinimum = projectedForce1 * projectedForce2;
            if (passedMinimum < 0. and (0.1 * fabs(projectedForce1) < fabs(projectedForce2)))
            {
                forceChange = (projectedForce1 - projectedForce2);
                stepSize = (projectedForce1 / forceChange) * stepSize;
            }
        }
    }
    if (parameters->optCGKnockOutMaxMove)
    {
        if (stepSize >= maxMove)
        {
            // knockout old search direction
            directionOld = objf->getPositions() * 0.0;
            forceOld = objf->getPositions() * 0.0;
        }
    }
    
//    return objf->isConverged();
    if(objf->isConverged()) return 1;
    return 0;
}


int ConjugateGradients::run(int maxIterations, double maxMove)
{
    int iterations = 0;
    while(!objf->isConverged() && iterations <= maxIterations)
    {
        step(maxMove);
        iterations++;
    }
//    return objf->isConverged();
    if(objf->isConverged()) return 1;
    return 0;
}


ConjugateGradients::~ConjugateGradients()
{
    return;
}

