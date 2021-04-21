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
    string posInFilename("pos.con");
    string posOutFilename("min.con");

    if (parameters->checkpoint) {
        FILE *pos;
        pos = fopen("pos_cp.con", "r");
        if (pos != NULL) {
            posInFilename = "pos_cp.con";
            log("[Minimization] Resuming from checkpoint\n");
        }else{
            log("[Minimization] No checkpoint files found\n");
        }
    }

    std::vector<std::string> returnFiles;
    returnFiles.push_back(posOutFilename);

    Matter *pos = new Matter(parameters);
    pos->con2matter(posInFilename);

    printf("\nBeginning minimization of %s\n", posInFilename.c_str());

    int status;

    bool converged;
    try {
        converged = pos->relax(false, parameters->writeMovies, 
                                    parameters->checkpoint, "minimization", "pos");
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

    printf("Saving result to %s\n", posOutFilename.c_str());
    pos->matter2con(posOutFilename);
    if (status != STATUS_POTENTIAL_FAILED) {
        printf("Final Energy: %f\n", pos->getPotentialEnergy());
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
        fprintf(fileResults, "%f potential_energy\n", pos->getPotentialEnergy());
    }
    fclose(fileResults);

    return returnFiles;
}
