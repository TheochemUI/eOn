//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include "SaddleSearchJob.h"
#include "Constants.h"
#include "ConjugateGradients.h"
#include "false_boinc.h"
#include "Potentials.h"

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

void SaddleSearchJob::run(int bundleNumber)
{
    char buff[STRING_SIZE];
    string reactant_passed("reactant_passed.con");
    string displacement_passed("displacement_passed");
    string mode_passed("mode_passed");

    if (bundleNumber < 0) {
        displacement_passed += ".con";
        mode_passed += ".dat";
    }else{
        snprintf(buff, STRING_SIZE, "_%i.con", bundleNumber);
        displacement_passed += buff;
        snprintf(buff, STRING_SIZE, "_%i.dat", bundleNumber);
        mode_passed += buff;
    }

    initial = new Matter(parameters);
    displacement = new Matter(parameters);
    saddle = new Matter(parameters);

    initial->con2matter(reactant_passed);
/*
    if (parameters->processSearchMinimizeFirst) {
        printf("Minimizing initial structure\n");
        int fi = Potential::fcalls;
        ConjugateGradients cgMin(initial, parameters);
        cgMin.fullRelax();
        fCallsMin += Potential::fcalls - fi;
    }
*/
    if (!parameters->saddleDisplaceType) {
        saddle->con2matter(displacement_passed);
    }
    else {
        *saddle = *initial;
    }
    saddlePoint = new SaddlePoint();
    saddlePoint->initialize(initial, saddle, parameters);
    if (!parameters->saddleDisplaceType) {
        saddlePoint->loadMode(mode_passed);
    }
    else {
        saddlePoint->displaceAndSetMode(saddle);
    }    

    int status = doSaddleSearch();

    printEndState(status);
    saveData(status, bundleNumber);

    delete saddlePoint;
    delete initial;
    delete displacement;
    delete saddle; 
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

void SaddleSearchJob::saveData(int status, int bundleNumber){
    FILE *fileResults, *fileSaddle, *fileMode;

    char filename[STRING_SIZE];

    if (bundleNumber != -1) {
        snprintf(filename, STRING_SIZE, "results_%i.dat", bundleNumber);
    }else{
        strncpy(filename, "results.dat", STRING_SIZE);
    }
	fileResults = fopen(filename, "wb");
	///XXX: min_fcalls isn't quite right it should get them from
	//      the minimizer. But right now the minimizers are in
	//      the SaddlePoint object. They will be taken out eventually.
    
    fprintf(fileResults, "%d termination_reason\n", status);
    fprintf(fileResults, "%ld random_seed\n", parameters->randomSeed);
    fprintf(fileResults, "%ld potential_tyep\n", parameters->potential);
    fprintf(fileResults, "%d total_force_calls\n", Potential::fcalls);
//    fprintf(fileResults, "%ld force_calls_minimization\n", saddlePoint->forceCallsMinimization + fCallsMin);
    fprintf(fileResults, "%d force_calls_saddle\n", fCallsSaddle);
    fprintf(fileResults, "%f potential_energy_saddle\n", saddle->getPotentialEnergy());
//    fprintf(fileResults, "%f potential_energy_reactant\n", min1->getPotentialEnergy());
//    fprintf(fileResults, "%f potential_energy_product\n", min2->getPotentialEnergy());
//    fprintf(fileResults, "%f barrier_reactant_to_product\n", barriersValues[0]);
//    fprintf(fileResults, "%f barrier_product_to_reactant\n", barriersValues[1]);
//    fprintf(fileResults, "%f displacement_saddle_distance\n",
//            displacement->perAtomNorm(*saddle));
//    fprintf(fileResults, "%d force_calls_prefactors\n", fCallsPrefactors);
//    fprintf(fileResults, "%.4e prefactor_reactant_to_product\n", prefactorsValues[0]);
//    fprintf(fileResults, "%.4e prefactor_product_to_reactant\n", prefactorsValues[1]);
	fclose(fileResults);

//    if (bundleNumber != -1) {
//        snprintf(filename, STRING_SIZE, "reactant_%i.con", bundleNumber);
//    }else{
//        strncpy(filename, "reactant.con", STRING_SIZE);
//    }
//	fileReactant = fopen(filename, "wb");
//    min1->matter2con(fileReactant);
//
    if (bundleNumber != -1) {
        snprintf(filename, STRING_SIZE, "mode_%i.dat", bundleNumber);
    }else{
        strncpy(filename, "mode.dat", STRING_SIZE);
    }
    fileMode = fopen(filename, "wb");
    saddlePoint->saveMode(fileMode);
    fclose(fileMode);
//	fclose(fileReactant);

    if (bundleNumber != -1) {
        snprintf(filename, STRING_SIZE, "saddle_%i.con", bundleNumber);
    }else{
        strncpy(filename, "saddle.con", STRING_SIZE);
    }
	fileSaddle = fopen(filename, "wb");
    saddle->matter2con(fileSaddle);
	fclose(fileSaddle);

//    if (bundleNumber != -1) {
//        snprintf(filename, STRING_SIZE, "product_%i.con", bundleNumber);
//    }else{
//        strncpy(filename, "product.con", STRING_SIZE);
//    }
//	fileProduct = fopen(filename, "wb");
//	min2->matter2con(fileProduct);
//    fclose(fileProduct);

    return;
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

//    else if(status == SaddlePoint::STATUS_BAD_NOT_CONNECTED)
//        fprintf(stdout, "Minima, saddle is not connected to initial state.\n");

//    else if(status == SaddlePoint::STATUS_BAD_PREFACTOR)
//            fprintf(stdout, "Prefactors, not within window\n");

//    else if(status == SaddlePoint::STATUS_BAD_HIGH_BARRIER)
//        fprintf(stdout, "Energy barriers, not within window\n");

//    else if (status == SaddlePoint::STATUS_BAD_MINIMA)
//        fprintf(stdout, "Minima, from saddle did not converge\n");
    else
        fprintf(stdout, "Unknown status: %i!\n", status);

    return;
}

