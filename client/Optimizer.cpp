#include "Optimizer.h"
#include "ConjugateGradients.h"
#include "Quickmin.h"
#include "LBFGS.h"
#include "FIRE.h"
#include "SteepestDescent.h"

Optimizer *Optimizer::getOptimizer(ObjectiveFunction *objf, Parameters *parameters)
{
    Optimizer* mizer = nullptr;
    if (parameters->optMethod == "cg") {
        mizer = dynamic_cast<Optimizer*>(new ConjugateGradients(objf, parameters));
    }else if (parameters->optMethod == "qm") {
        mizer = dynamic_cast<Optimizer*>(new Quickmin(objf, parameters));
    }else if (parameters->optMethod == "lbfgs") {
        mizer = dynamic_cast<Optimizer*>(new LBFGS(objf, parameters));
    }else if (parameters->optMethod == "fire") {
        mizer = dynamic_cast<Optimizer*>(new FIRE(objf, parameters));
    }else if (parameters->optMethod == "sd") {
        mizer = dynamic_cast<Optimizer*>(new SteepestDescent(objf, parameters));
    }else{
        printf("Unknown optMethod: %s\n", parameters->optMethod.c_str());
        exit(1);
    }
    return mizer;
}
