
#include "Quickmin.h"
#include "HelperFunctions.h"
#include <cmath>

Quickmin::Quickmin(ObjectiveFunction *objfPassed, Parameters *parametersPassed)
{
    objf = objfPassed;
    parameters = parametersPassed;
    dt = parametersPassed->optTimeStep;
    velocity.resize(objf->degreesOfFreedom());
    velocity.setZero();
    iteration = 0;
}

Quickmin::~Quickmin()
{
    return;
}

bool Quickmin::step(double maxMove)
{
    VectorXd force = -objf->getGradient();
    if (parameters->optQMSteepestDecent) {
        velocity.setZero();
    }
    else {
        if (velocity.dot(force) < 0) {
            velocity.setZero();
        }
        else {
            VectorXd f_unit = force/force.norm();
            velocity = velocity.dot(f_unit) * f_unit;
        }
    }
    
    velocity += force * dt;
    VectorXd dr = helper_functions::maxAtomMotionAppliedV(velocity * dt, parameters->optMaxMove);
    objf->setPositions(objf->getPositions() + dr);  
    iteration++;
    return objf->isConverged();
}

bool Quickmin::run(int maxSteps, double maxMove)
{
    while(!objf->isConverged() && iteration < maxSteps) {
        step(maxMove);
    }
    return objf->isConverged();
}
