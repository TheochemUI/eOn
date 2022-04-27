
//Based on the SteepestDescent minimizer written in ASE.

#include "SteepestDescent.h"
#include "Log.h"
#include <cassert>
#include <cmath>

SteepestDescent::SteepestDescent(ObjectiveFunction *objfPassed, Parameters *parametersPassed)
{
    objf = objfPassed;
    parameters = parametersPassed;

    iteration = 0;
}

SteepestDescent::~SteepestDescent()
{
    return;
}

int SteepestDescent::step(double maxMove)
{
    VectorXd r = objf->getPositions();
    VectorXd f = -objf->getGradient();

    VectorXd dr;
    double alpha = parameters->optSDAlpha;
    if (parameters->optSDTwoPoint == true && iteration > 0) {
        VectorXd dx = r-rPrev;
        VectorXd dg = -f+fPrev;
        alpha = dx.dot(dx)/dx.dot(dg);
        if (alpha < 0) {
            alpha = parameters->optSDAlpha;
        }
        log_file("[SD] alpha: %.4e\n", alpha);
    }

    dr = alpha*f;
    dr = helper_functions::maxAtomMotionAppliedV(dr, maxMove);

    objf->setPositions(r+dr);

    rPrev = r;
    fPrev = f;

    iteration++;

//    return objf->isConverged();
    if(objf->isConverged()) return 1;
    return 0;
}


int SteepestDescent::run(int maxSteps, double maxMove)
{
    while(!objf->isConverged() && iteration < maxSteps) {
        step(maxMove);
    }
//    return objf->isConverged();
    if(objf->isConverged()) return 1;
    return 0;
}

int SteepestDescent::step(const double maxMove,
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
