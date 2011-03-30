//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include "NudgedElasticBandJob.h"
#include "Constants.h"
#include "ConjugateGradients.h"
#include "false_boinc.h"
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
    string reactant_passed("reactant_passed.con");
    string displacement_passed("displacement_passed.con");
    string mode_passed("mode_passed.dat");

    initial = new Matter(parameters);
    displacement = new Matter(parameters);
    saddle = new Matter(parameters);

    initial->con2matter(reactant_passed);

    if (parameters->saddleDisplaceType == SaddlePoint::DISP_LOAD) {
        // displacement was passed from the server        
        saddle->con2matter(displacement_passed);
    }
    else {
        // displacement and mode will be made on the client
        // in saddlePoint->initialize(...) 
        *saddle = *initial;
    }
    saddlePoint = new SaddlePoint();
    saddlePoint->initialize(initial, saddle, parameters);
    if (parameters->saddleDisplaceType == SaddlePoint::DISP_LOAD) {
        // mode was passed from the server        
        saddlePoint->loadMode(mode_passed);
    }

    int status = doSaddleSearch();

    printEndState(status);
    saveData(status);

    delete saddlePoint;
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
    status = saddlePoint->locate();
    fCallsSaddle += Potential::fcalls - f1;

    if (status == SaddlePoint::STATUS_INIT) {
        status = SaddlePoint::STATUS_GOOD;
    }

    return status;
}

void SaddleSearchJob::saveData(int status){
    FILE *fileResults, *fileSaddle, *fileMode;

    std::string resultsFilename("results.dat");
    returnFiles.push_back(resultsFilename);
    fileResults = fopen(resultsFilename.c_str(), "wb");
    ///XXX: min_fcalls isn't quite right it should get them from
    //      the minimizer. But right now the minimizers are in
    //      the SaddlePoint object. They will be taken out eventually.
    
    fprintf(fileResults, "%d termination_reason\n", status);
    fprintf(fileResults, "%ld random_seed\n", parameters->randomSeed);
    fprintf(fileResults, "%s potential_tyep\n", parameters->potential.c_str());
    fprintf(fileResults, "%d total_force_calls\n", Potential::fcalls);
    fprintf(fileResults, "%d force_calls_saddle\n", fCallsSaddle);
    fprintf(fileResults, "%f potential_energy_saddle\n", saddle->getPotentialEnergy());
    fclose(fileResults);

    std::string modeFilename("mode.dat");
    returnFiles.push_back(modeFilename);
    fileMode = fopen(modeFilename.c_str(), "wb");
    saddlePoint->saveMode(fileMode);
    fclose(fileMode);

    std::string saddleFilename("saddle.con");
    returnFiles.push_back(saddleFilename);
    fileSaddle = fopen(saddleFilename.c_str(), "wb");
    saddle->matter2con(fileSaddle);
    fclose(fileSaddle);
}

void SaddleSearchJob::printEndState(int status) {
    fprintf(stdout, "Final state: ");
    if(status == SaddlePoint::STATUS_GOOD)
        fprintf(stdout, "Successful.\n");

    else if(status == SaddlePoint::STATUS_BAD_NO_CONVEX)
        fprintf(stdout, "Initial displacement, not able to reach convex region.\n");
   
    else if(status == SaddlePoint::STATUS_BAD_HIGH_ENERGY)
        fprintf(stdout, "Saddle search, barrier too high.\n");
        
    else if(status == SaddlePoint::STATUS_BAD_MAX_CONCAVE_ITERATIONS) 
        fprintf(stdout, "Saddle search, too many iterations in concave region.\n");

    else if(status == SaddlePoint::STATUS_BAD_MAX_ITERATIONS)
        fprintf(stdout, "Saddle search, too many iterations in saddle point search.\n");
    else
        fprintf(stdout, "Unknown status: %i!\n", status);

    return;
}

