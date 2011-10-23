//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include "HessianJob.h"
#include "Matter.h"
#include "Hessian.h"
#include "Potential.h"

HessianJob::HessianJob(Parameters *params)
{
    parameters = params;
}

HessianJob::~HessianJob()
{
}

std::vector<std::string> HessianJob::run(void)
{
    string reactant_passed("reactant_passed.con");
    string saddle_passed("saddle_passed.con");
    string product_passed("product_passed.con");

    std::vector<std::string> returnFiles;

    Matter *reactant = new Matter(parameters);
    Matter *saddle = new Matter(parameters);
    Matter *product = new Matter(parameters);

    reactant->con2matter(reactant_passed);
    saddle->con2matter(saddle_passed);
    product->con2matter(product_passed);

// GH: need to fix this

//    Hessian hessian(reactant, saddle, product, parameters);
//    VectorXd modes = hessian.getModes(parameters->hessianType);
    //XXX: to fix build for now...
    VectorXd modes;

    bool failed = modes.size()==0;

    FILE *fileResults;
    FILE *fileMode;

    std::string results_file("results.dat");
    std::string mode_file("mode.dat");

    returnFiles.push_back(results_file);
    returnFiles.push_back(mode_file);

    fileResults = fopen(results_file.c_str(), "wb");
    fileMode = fopen(mode_file.c_str(), "wb");

    fprintf(fileResults, "%s good\n", failed ? "false" : "true");
    fprintf(fileResults, "%d force_calls\n", Potential::fcalls);
//XXX: to fix build
//    fprintf(fileResults, "%d hessian_size\n", 
//            (int)hessian.getHessian(parameters->hessianType).rows());
    if(!failed)
    {
        for(int i=0; i<modes.size(); i++)
        {
            fprintf(fileMode, "%f\n", modes[i]);
        }
    }

    delete reactant;
    delete product;
    delete saddle;

    return returnFiles;
}
 
