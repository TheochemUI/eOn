
#include "SaddleSearchJob.h"
#include "EpiCenters.h"
#include "Log.h"
#include "Potential.h"

#include <stdio.h>
#include <string>

using namespace std;

SaddleSearchJob::SaddleSearchJob(Parameters *params)
{
    parameters = params;
    fCallsSaddle = 0;
}

SaddleSearchJob::~SaddleSearchJob()
{}

std::vector<std::string> SaddleSearchJob::run(void)
{
    string reactantFilename("pos.con");
    string displacementFilename("displacement.con");
    string modeFilename("direction.dat");

    if (parameters->checkpoint) {
        FILE *disp, *mode;
        disp = fopen("displacement_cp.con", "r");
        mode = fopen("mode_cp.dat", "r");
        if (disp != NULL && mode != NULL) {
            displacementFilename = "displacement_cp.con";
            modeFilename = "mode_cp.dat";
            log("Resuming from checkpoint\n");
        }else{
            log("No checkpoint files found\n");
        }
    }

    initial = new Matter(parameters);
    displacement = new Matter(parameters);
    saddle = new Matter(parameters);

    initial->con2matter(reactantFilename);

    if (parameters->saddleDisplaceType == EpiCenters::DISP_LOAD) {
        // displacement was passed from the server
        saddle->con2matter(displacementFilename);
    }
    else {
        // displacement and mode will be made on the client
        // in saddleSearch->initialize(...)
        *saddle = *initial;
    }
    AtomMatrix mode;
    if (parameters->saddleDisplaceType == EpiCenters::DISP_LOAD) {
        // mode was passed from the server
        mode = helper_functions::loadMode(modeFilename, initial->numberOfAtoms());
    }

    saddleSearch = new MinModeSaddleSearch(saddle, mode, initial->getPotentialEnergy(), parameters);

    int status;
    status = doSaddleSearch();
    printEndState(status);
    saveData(status);

    delete saddleSearch;
    delete initial;
    delete displacement;
    delete saddle; 

    return returnFiles;
}

int SaddleSearchJob::doSaddleSearch()
{
    Matter matterTemp(parameters);
    long status;
    int f1;
    f1 = Potential::fcalls;
    try {
        status = saddleSearch->run();
    }catch (int e) {
        if (e == 100) {
            status = MinModeSaddleSearch::STATUS_POTENTIAL_FAILED; 
        }else{
            printf("unknown exception: %i\n", e);
            throw e;
        }
    }

    fCallsSaddle += Potential::fcalls - f1;

    return status;
}

void SaddleSearchJob::saveData(int status){
    FILE *fileResults, *fileSaddle, *fileMode;

    std::string resultsFilename("results.dat");
    returnFiles.push_back(resultsFilename);
    fileResults = fopen(resultsFilename.c_str(), "wb");
    ///XXX: min_fcalls isn't quite right it should get them from
    //      the minimizer. But right now the minimizers are in
    //      the SaddleSearch object. They will be taken out eventually.

    fprintf(fileResults, "%d termination_reason\n", status);
    fprintf(fileResults, "saddle_search job_type\n");
    fprintf(fileResults, "%ld random_seed\n", parameters->randomSeed);
    fprintf(fileResults, "%s potential_type\n", parameters->potential.c_str());
    fprintf(fileResults, "%d total_force_calls\n", Potential::fcalls);
    fprintf(fileResults, "%d force_calls_saddle\n", fCallsSaddle);
    fprintf(fileResults, "%i iterations\n", saddleSearch->iteration);
    if (status != MinModeSaddleSearch::STATUS_POTENTIAL_FAILED) {
        fprintf(fileResults, "%f potential_energy_saddle\n", saddle->getPotentialEnergy());
        fprintf(fileResults, "%f final_eigenvalue\n", saddleSearch->getEigenvalue());
    }
    fprintf(fileResults, "%f potential_energy_reactant\n", initial->getPotentialEnergy());
    fclose(fileResults);

    std::string modeFilename("mode.dat");
    returnFiles.push_back(modeFilename);
    fileMode = fopen(modeFilename.c_str(), "wb");
    helper_functions::saveMode(fileMode, saddle, saddleSearch->getEigenvector());
    fclose(fileMode);

    std::string saddleFilename("saddle.con");
    returnFiles.push_back(saddleFilename);
    fileSaddle = fopen(saddleFilename.c_str(), "wb");
    saddle->matter2con(fileSaddle);
    fclose(fileSaddle);
}

void SaddleSearchJob::printEndState(int status)
{
    log("[Saddle Search] Final status: ");

    if(status == MinModeSaddleSearch::STATUS_GOOD)
        log("Success\n");

    else if(status == MinModeSaddleSearch::STATUS_BAD_NO_CONVEX)
        log("Initial displacement unable to reach convex region\n");

    else if(status == MinModeSaddleSearch::STATUS_BAD_HIGH_ENERGY)
        log("Barrier too high\n");

    else if(status == MinModeSaddleSearch::STATUS_BAD_MAX_CONCAVE_ITERATIONS)
        log("Too many iterations in concave region\n");

    else if(status == MinModeSaddleSearch::STATUS_BAD_MAX_ITERATIONS)
        log("Too many iterations\n");

    else if(status == MinModeSaddleSearch::STATUS_BAD_NOT_CONNECTED)
        log("Saddle is not connected to initial state\n");

    else if(status == MinModeSaddleSearch::STATUS_BAD_PREFACTOR)
        log("Prefactors not within window\n");

    else if(status == MinModeSaddleSearch::STATUS_FAILED_PREFACTOR)
        log("Hessian calculation failed\n");

    else if(status == MinModeSaddleSearch::STATUS_BAD_HIGH_BARRIER)
        log("Energy barrier not within window\n");

    else if (status == MinModeSaddleSearch::STATUS_BAD_MINIMA)
        log("Minimizations from saddle did not converge\n");

    else if(status == MinModeSaddleSearch::STATUS_NONNEGATIVE_ABORT)
        log("Nonnegative initial mode, aborting\n");

    else if(status == MinModeSaddleSearch::STATUS_NEGATIVE_BARRIER)
        log("Negative barrier detected\n");

    else if(status == MinModeSaddleSearch::STATUS_BAD_MD_TRAJECTORY_TOO_SHORT)
        log("No reaction found during MD trajectory\n");

    else if(status == MinModeSaddleSearch::STATUS_BAD_NO_NEGATIVE_MODE_AT_SADDLE)
        log("Converged to stationary point with zero negative modes\n");

    else if(status == MinModeSaddleSearch::STATUS_BAD_NO_BARRIER)
        log("No forward barrier was found along minimized band\n");

    else if(status == MinModeSaddleSearch::STATUS_ZEROMODE_ABORT)
        log("Zero mode abort.\n");

    else if(status == MinModeSaddleSearch::STATUS_OPTIMIZER_ERROR)
        log("Optimizer error.\n");

    else
        log("Unknown status: %i!\n", status);

    return;
}
