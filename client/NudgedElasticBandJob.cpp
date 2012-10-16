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
#include "Log.h"

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
    long status;
    int f1;

    string reactantFilename = helper_functions::getRelevantFile("reactant.con");
    string productFilename = helper_functions::getRelevantFile("product.con");

    Matter *initial = new Matter(parameters);
    Matter *final = new Matter(parameters);

    initial->con2matter(reactantFilename);
    final->con2matter(productFilename);

    NudgedElasticBand *neb = new NudgedElasticBand(initial, final, parameters);

    f1 = Potential::fcalls;
    status = neb->compute();
    fCallsNEB += Potential::fcalls - f1;

    if (status == NudgedElasticBand::STATUS_INIT) {
        status = NudgedElasticBand::STATUS_GOOD;
    }

    printEndState(status);
    saveData(status, neb);

    delete neb;
    delete initial;
    delete final;

    return returnFiles;
}

void NudgedElasticBandJob::saveData(int status, NudgedElasticBand *neb)
{
    FILE *fileResults, *fileNEB;

    std::string resultsFilename("results.dat");
    returnFiles.push_back(resultsFilename);
    fileResults = fopen(resultsFilename.c_str(), "wb");
    
    fprintf(fileResults, "%d termination_reason\n", status);
    fprintf(fileResults, "%s potential_type\n", parameters->potential.c_str());
    fprintf(fileResults, "%d total_force_calls\n", Potential::fcalls);
    fprintf(fileResults, "%d force_calls_neb\n", fCallsNEB);
    fprintf(fileResults, "%f energy_reference\n", neb->image[0]->getPotentialEnergy());
    fprintf(fileResults, "%li number_of_images\n", neb->images);
    for(long i=0; i<=neb->images+1; i++) {
        fprintf(fileResults, "%f image%li_energy\n", neb->image[i]->getPotentialEnergy()-neb->image[0]->getPotentialEnergy(), i);
        fprintf(fileResults, "%f image%li_force\n", neb->image[i]->getForces().norm(), i);
        fprintf(fileResults, "%f image%li_projected_force\n", neb->projectedForce[i]->norm(), i);
    }
    fprintf(fileResults, "%li number_of_extrema\n", neb->numExtrema);
    for(long i=0; i<neb->numExtrema; i++) {
        fprintf(fileResults, "%f extremum%li_position\n", neb->extremumPosition[i], i);
        fprintf(fileResults, "%f extremum%li_energy\n", neb->extremumEnergy[i], i);
    }

    fclose(fileResults);

    std::string nebFilename("neb.con");
    returnFiles.push_back(nebFilename);
    fileNEB = fopen(nebFilename.c_str(), "wb");
    for(long i=0; i<=neb->images+1; i++) {
        neb->image[i]->matter2con(fileNEB);
    }
    fclose(fileNEB);
}

void NudgedElasticBandJob::printEndState(int status)
{
    log("\nFinal state: ");
    if(status == NudgedElasticBand::STATUS_GOOD)
        log("Nudged elastic band, successful.\n");
    else if(status == NudgedElasticBand::STATUS_BAD_MAX_ITERATIONS)
        log("Nudged elastic band, too many iterations.\n");
    else
        log("Unknown status: %i!\n", status);
    return;
}

