#include "MonteCarloJob.h"
#include "MonteCarlo.h"
#include "Log.h"
#include "Matter.h"
#include "HelperFunctions.h"

MonteCarloJob::MonteCarloJob(Parameters *params)
{
    parameters = params;
}

MonteCarloJob::~MonteCarloJob(){ }

std::vector<std::string> MonteCarloJob::run(void)
{
    string posInFilename("pos.con");
    string posOutFilename("out.con");

    if (parameters->checkpoint) {
        FILE *pos;
        pos = fopen("pos_cp.con", "r");
        if (pos != NULL) {
            posInFilename = "pos_cp.con";
            log("Resuming from checkpoint\n");
        }else{
            log("No checkpoint files found\n");
        }
    }

    std::vector<std::string> returnFiles;
    returnFiles.push_back(posOutFilename);

    Matter *matter = new Matter(parameters);
    matter->con2matter(posInFilename);

    //code will go
    MonteCarlo mc = MonteCarlo(matter, parameters);
    mc.run(parameters->monteCarloSteps, parameters->temperature, 
            parameters->monteCarloStepSize);


    //FILE *fileResults;

    std::string resultsFilename("results.dat");
    returnFiles.push_back(resultsFilename);
    //fileResults = fopen(resultsFilename.c_str(), "wb");

    //fprintf(fileResults, "%d termination_reason\n", status);
    //fprintf(fileResults, "minimization job_type\n");
    //fprintf(fileResults, "%s potential_type\n", parameters->potential.c_str());
    //fprintf(fileResults, "%d total_force_calls\n", Potential::fcallsTotal);
    //if (status != STATUS_POTENTIAL_FAILED) {
    //    fprintf(fileResults, "%f potential_energy\n", pos->getPotentialEnergy());
    //}
    //fclose(fileResults);

    return returnFiles;
}
