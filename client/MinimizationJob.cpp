//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include "MinimizationJob.h"
#include "Minimizer.h"
#include "ConjugateGradients.h"
#include "Quickmin.h"
#include "QuickminBox.h"
#include "Matter.h"
#include "Constants.h"

MinimizationJob::MinimizationJob(Parameters *params)
{
    parameters = params;
}

MinimizationJob::~MinimizationJob(){ }

std::vector<std::string> MinimizationJob::run(void)
{
    string reactant_passed("reactant_passed.con");
    string reactant_output("reactant.con");

    std::vector<std::string> returnFiles;
    returnFiles.push_back(reactant_output);

    Matter *reactant = new Matter(parameters);
    reactant->con2matter(reactant_passed);

    printf("\nBeginning minimization of %s\n", reactant_passed.c_str());

    int status;

    Minimizer* mizer=NULL;
    try {
        if(parameters->optMethod == "cg")
        {
            mizer = new ConjugateGradients(reactant, parameters);
        }
        else if(parameters->optMethod == "qm")
        {
            mizer = new Quickmin(reactant, parameters);
        }
        else if(parameters->optMethod == "box")
        {
            mizer = new QuickminBox(reactant, parameters);
        }else{
            printf("Unknown optMethod: %s\n", parameters->optMethod.c_str());
        }
    }catch (int e) {
        if (e == 100) {
            status = MinimizationJob::STATUS_POTENTIAL_FAILED;
        }else{
            printf("unknown exception: %i\n", e);
            throw e;
        }
    }

    mizer->setOutput(1);
    mizer->fullRelax();
    if (mizer->isItConverged(parameters->optConvergedForce)) {
        printf("Minimization converged within tolerence\n");
        status = MinimizationJob::STATUS_GOOD;
    }else{
        printf("Minimization did not converge to tolerence!\n"
               "Maybe try to increase maximum_iterations?\n");
        status = MinimizationJob::STATUS_MAX_ITERATIONS;
    }

    printf("Saving result to %s\n", reactant_output.c_str());
    reactant->matter2con(reactant_output);
    if (status != MinimizationJob::STATUS_POTENTIAL_FAILED) {
        printf("Final Energy: %f\n", reactant->getPotentialEnergy());
    }

    FILE *fileResults;

    std::string resultsFilename("results.dat");
    returnFiles.push_back(resultsFilename);
    fileResults = fopen(resultsFilename.c_str(), "wb");

    fprintf(fileResults, "%d termination_reason\n", status);
    fprintf(fileResults, "minimization job_type\n");
    fprintf(fileResults, "%s potential_type\n", parameters->potential.c_str());
    fprintf(fileResults, "%d total_force_calls\n", Potential::fcalls);
    if (status != MinimizationJob::STATUS_POTENTIAL_FAILED) {
        fprintf(fileResults, "%f potential_energy\n", reactant->getPotentialEnergy());
    }
    fclose(fileResults);

    return returnFiles;
}
