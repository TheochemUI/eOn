#include "Optimizer.h"
#include "ConjugateGradients.h"
#include "Quickmin.h"
#include "LBFGS.h"
#include "FIRE.h"
#include "SteepestDescent.h"

Optimizer *Optimizer::getOptimizer(ObjectiveFunction *objf, Parameters *parameters)
{
    Optimizer* optimizer = NULL;
    if (parameters->optMethod == "cg") {
        optimizer = new ConjugateGradients(objf, parameters);
    } else if (parameters->optMethod == "qm") {
        optimizer = new Quickmin(objf, parameters);
    } else if (parameters->optMethod == "lbfgs") {
        optimizer = new LBFGS(objf, parameters);
    } else if (parameters->optMethod == "fire") {
        optimizer = new FIRE(objf, parameters);
    } else if (parameters->optMethod == "sd") {
        optimizer = new SteepestDescent(objf, parameters);
    } else {
        printf("Unknown optMethod: %s\n", parameters->optMethod.c_str());
        exit(1);
    }
    return optimizer;
}
