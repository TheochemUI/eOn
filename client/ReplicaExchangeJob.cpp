#include <cstdlib>

#include "Matter.h"
#include "Dynamics.h"
#include "BondBoost.h"
#include "ReplicaExchangeJob.h"
#include "ConjugateGradients.h"

ReplicaExchangeJob::ReplicaExchangeJob(Parameters *params)
{
    parameters = params;
}

ReplicaExchangeJob::~ReplicaExchangeJob(){ }

void ReplicaExchangeJob::run(int bundleNumber)
{
    char buff[STRING_SIZE];
    string reactant_passed("reactant_passed");

    if (bundleNumber < 0) {
        reactant_passed += ".con";
    }else{
        snprintf(buff, STRING_SIZE, "_%i.con", bundleNumber);
        reactant_passed += buff;
    }

    reactant = new Matter(parameters);

    reactant->con2matter(reactant_passed);

    printf("Now running Parralel Replica Dynamics\n");

    delete reactant;
}

