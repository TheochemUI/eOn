#include <stdio.h>
#include <string>

#include "Dynamics.h"
#include "DynamicsJob.h"
#include "Parameters.h"
#include "Potential.h"
#include "HelperFunctions.h"

using namespace helper_functions;

std::vector<std::string> DynamicsJob::run(void)
{
    Matter *R = new Matter(params);
    Matter *F = new Matter(params);
    R->con2matter("pos.con");
    *F = *R;

    Dynamics *d = new Dynamics(R, params.get());
    d->run();

    *F = *R;
    FILE *fileProduct;
    std::string productFilename("final.con");
    returnFiles.push_back(productFilename);

    fileProduct = fopen(productFilename.c_str(), "wb");
    F->matter2con(fileProduct);
    fclose(fileProduct);

    delete R;
    delete F;
    delete d;

    std::vector<std::string> returnFiles;
    return returnFiles;
}

