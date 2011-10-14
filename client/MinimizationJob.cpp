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
#include "Optimizer.h"
#include "Log.h"
#include "Matter.h"
#include "HelperFunctions.h"

MinimizationJob::MinimizationJob(Parameters *params)
{
    parameters = params;
    fcalls = Potential::fcalls;
}

MinimizationJob::~MinimizationJob(){ }

std::vector<std::string> MinimizationJob::run(void)
{
//    string reactant_passed("reactant_passed.con");
    string reactant_passed = helper_functions::getRelevantFile(parameters->conFilename);
    string reactant_output("reactant.con");

    if (parameters->checkpoint) {
        FILE *react;
        react = fopen("reactant_checkpoint.con", "r");
        if (react != NULL) {
            reactant_passed = "reactant_checkpoint.con";
            log("Resuming from checkpoint\n");
        }else{
            log("No checkpoint files found\n");
        }
    }

    std::vector<std::string> returnFiles;
    returnFiles.push_back(reactant_output);

    Matter *reactant = new Matter(parameters);
    reactant->con2matter(reactant_passed);

    printf("\nBeginning minimization of %s\n", reactant_passed.c_str());


    int status;

    bool converged;
    try {
        converged = reactant->relax(false, parameters->writeMovies, 
                                    parameters->checkpoint, "min", "reactant");
        if (converged) {
            status = STATUS_GOOD;
            printf("Minimization converged within tolerence\n");
        }else{
            status = STATUS_MAX_ITERATIONS;
            printf("Minimization did not converge to tolerence!\n"
                   "Maybe try to increase max_iterations?\n");
        }
    }catch (int e) {
        if (e == 100) {
            status = STATUS_POTENTIAL_FAILED;
        }else{
            throw e;
        }
    }

    printf("Saving result to %s\n", reactant_output.c_str());
    reactant->matter2con(reactant_output);
    if (status != STATUS_POTENTIAL_FAILED) {
        printf("Final Energy: %f\n", reactant->getPotentialEnergy());
    }

    FILE *fileResults;

    std::string resultsFilename("results.dat");
    returnFiles.push_back(resultsFilename);
    fileResults = fopen(resultsFilename.c_str(), "wb");

    fprintf(fileResults, "%d termination_reason\n", status);
    fprintf(fileResults, "minimization job_type\n");
    fprintf(fileResults, "%s potential_type\n", parameters->potential.c_str());
    fprintf(fileResults, "%d total_force_calls\n", Potential::fcallsTotal);
    if (status != STATUS_POTENTIAL_FAILED) {
        fprintf(fileResults, "%f potential_energy\n", reactant->getPotentialEnergy());
    }
    fclose(fileResults);

    return returnFiles;
}
