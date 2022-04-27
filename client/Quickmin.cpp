
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

int Quickmin::step(double maxMove)
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
//    return objf->isConverged();
    if(objf->isConverged()) return 1;
    return 0;
}

int Quickmin::run(int maxSteps, double maxMove)
{
    while(!objf->isConverged() && iteration < maxSteps) {
        step(maxMove);
    }
//    return objf->isConverged();
    if(objf->isConverged()) return 1;
    return 0;
}

int Quickmin::step(const double maxMove,
                             const std::vector<Matter> ppoints,
                             const double max_dist,
                             bool& isWithin){
    int stepval = step(maxMove);
    size_t nfree = ppoints.front().numberOfFreeAtoms();
    VectorXd cpath = this->objf->getPositions();
    // We also need to calculate how many objects are present in cpath
    // TODO: This assumes only free atoms are in getPositions
    size_t nimg = cpath.size() / 3 * nfree;
    // We need to construct a false Matter object to ensure distances respect PBCs
    Matter tmpm = ppoints.front();
    for (size_t idx{1}; idx < nimg-1; idx++){
        tmpm.setPositionsFreeV(
            cpath.segment(3 * nfree * (idx-1), 3 * nfree)
            );
        std::vector<double> distances;
        std::transform(ppoints.begin(),
                       ppoints.end(),
                       std::back_inserter(distances),
                       [&](Matter mat)->double{
                           return mat.distanceTo(tmpm);
                       });
        distances.erase(std::remove(distances.begin(), distances.end(), 0), distances.end());
        isWithin = std::any_of(distances.begin(),
                          distances.end(),
                          [&](const double dist)->bool{
                              return (std::abs(dist) < max_dist); });
        if (!isWithin){
            std::cout<<"\nEARLY STOPPING in OPTIMIZER\n";
            break;
        }
    }
    return stepval;
}
