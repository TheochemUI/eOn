#include "Optimizer.h"
#include "ConjugateGradients.h"
#include "Quickmin.h"
#include "LBFGS.h"
#include "FIRE.h"
#include "SteepestDescent.h"

Optimizer *Optimizer::getOptimizer(ObjectiveFunction *objf, Parameters *parameters)
{
    Optimizer* mizer=NULL;
    if (parameters->optMethod == "cg") {
        mizer = new ConjugateGradients(objf, parameters);
    }else if (parameters->optMethod == "qm") {
        mizer = new Quickmin(objf, parameters);
    }else if (parameters->optMethod == "lbfgs") {
        mizer = new LBFGS(objf, parameters);
    }else if (parameters->optMethod == "fire") {
        mizer = new FIRE(objf, parameters);
    }else if (parameters->optMethod == "sd") {
        mizer = new SteepestDescent(objf, parameters);
    }else{
        printf("Unknown optMethod: %s\n", parameters->optMethod.c_str());
        exit(1);
    }
    return mizer;
}
