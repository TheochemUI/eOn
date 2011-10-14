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
#include "ConjugateGradients.h"
#include "false_boinc.h"
#include "Potential.h"

#include <stdio.h>
#include <string>

using namespace std;

NudgedElasticBandJob::NudgedElasticBandJob(Parameters *parametersPassed)
{
    parameters = parametersPassed;
    fCallsNEB = 0;
}

NudgedElasticBandJob::~NudgedElasticBandJob()
{}

std::vector<std::string> NudgedElasticBandJob::run(void)
{
    string reactant_passed("reactant_passed.con");
    string product_passed("product_passed.con");

    initial = new Matter(parameters);
    final = new Matter(parameters);

    initial->con2matter(reactant_passed);
    final->con2matter(product_passed);

    neb = new NudgedElasticBand(initial, final, parameters);

    int status = findMinimumEnergyPath();

    printEndState(status);
    saveData(status);

    delete neb;
    delete initial;
    delete final;

    return returnFiles;
}

int NudgedElasticBandJob::findMinimumEnergyPath()
{
    Matter matterTemp(parameters);
    long status;
    int f1;
    f1 = Potential::fcalls;
    status = neb->compute();
    fCallsNEB += Potential::fcalls - f1;

    if (status == NudgedElasticBand::STATUS_INIT) {
        status = NudgedElasticBand::STATUS_GOOD;
    }

    return status;
}

void NudgedElasticBandJob::saveData(int status){
    FILE *fileResults, *fileNEB, *fileMode;

    std::string resultsFilename("results.dat");
    returnFiles.push_back(resultsFilename);
    fileResults = fopen(resultsFilename.c_str(), "wb");
    ///XXX: min_fcalls isn't quite right it should get them from
    //      the minimizer. But right now the minimizers are in
    //      the NudgedElasticBand object. They will be taken out eventually.
    
    fprintf(fileResults, "%d termination_reason\n", status);
    fprintf(fileResults, "%s potential_type\n", parameters->potential.c_str());
    fprintf(fileResults, "%d total_force_calls\n", Potential::fcalls);
    fprintf(fileResults, "%d force_calls_neb\n", fCallsNEB);
//    fprintf(fileResults, "%f potential_energy_saddle\n", neb->getPotentialEnergy());
    fclose(fileResults);

/*
    std::string saddleFilename("saddle.con");
    returnFiles.push_back(saddleFilename);
    fileSaddle = fopen(saddleFilename.c_str(), "wb");
    saddle->matter2con(fileSaddle);
    fclose(fileSaddle);
*/
}

void NudgedElasticBandJob::printEndState(int status) {
/*    fprintf(stdout, "Final state: ");
    if(status == NudgedElasticBand::STATUS_GOOD)
        fprintf(stdout, "Successful.\n");

    else if(status == NudgedElasticBand::STATUS_BAD_NO_CONVEX)
        fprintf(stdout, "Initial displacement, not able to reach convex region.\n");
   
    else if(status == NudgedElasticBand::STATUS_BAD_HIGH_ENERGY)
        fprintf(stdout, "Saddle search, barrier too high.\n");
        
    else if(status == NudgedElasticBand::STATUS_BAD_MAX_CONCAVE_ITERATIONS) 
        fprintf(stdout, "Saddle search, too many iterations in concave region.\n");

    else if(status == NudgedElasticBand::STATUS_BAD_MAX_ITERATIONS)
        fprintf(stdout, "Saddle search, too many iterations in saddle point search.\n");
    else
        fprintf(stdout, "Unknown status: %i!\n", status);
*/
    return;
}

