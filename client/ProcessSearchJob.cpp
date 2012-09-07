//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include "ProcessSearchJob.h"
#include "EpiCenters.h"
#include "Optimizer.h"
#include "false_boinc.h"
#include "Potential.h"
#include "Log.h"
#include "Prefactor.h"
#include "DynamicsSaddleSearch.h"

#include <stdio.h>
#include <string>
//#include <cassert>

ProcessSearchJob::ProcessSearchJob (Parameters *params)
{
    parameters = params;
    fCallsSaddle = fCallsPrefactors = fCallsMin = 0;
}

ProcessSearchJob::~ProcessSearchJob()
{}

std::vector<std::string> ProcessSearchJob::run(void)
{
    string reactantFilename("pos_in.con");
    string displacementFilename("displacement.con");
    string modeFilename("mode_in.dat");

    initial = new Matter(parameters);
    if (parameters->saddleMethod == "min_mode") {
        displacement = new Matter(parameters);
    }else{
        displacement = NULL;
    }
    saddle = new Matter(parameters);
    min1 = new Matter(parameters);
    min2 = new Matter(parameters);

    if (!initial->con2matter(reactantFilename)) {
        printf("Stop\n");
        exit(1);
    }

    if (parameters->processSearchMinimizeFirst) {
        log("Minimizing initial structure\n");
        int fi = Potential::fcalls;
        initial->relax();
        fCallsMin += Potential::fcalls - fi;
    }

    barriersValues[0] = barriersValues[1] = 0;
    prefactorsValues[0] = prefactorsValues[1] = 0;

    if (parameters->saddleMethod == "min_mode") {
        if (parameters->saddleDisplaceType == EpiCenters::DISP_LOAD) {
            // displacement was passed from the server
            if(!saddle->con2matter(displacementFilename)) {
                printf("Stop\n");
                exit(1);
            }
            *min1 = *min2 = *initial;
        }else{
            // displacement and mode will be made on the client
            // in SaddleSearch->initialize(...)
            *saddle = *min1 = *min2 = *initial;
        }
    }else{
        *saddle = *min1 = *min2 = *initial;
    }

    AtomMatrix mode;

    if (parameters->saddleMethod == "min_mode") {
        if (parameters->saddleDisplaceType == EpiCenters::DISP_LOAD) {
            // mode was passed from the server
            mode = helper_functions::loadMode(modeFilename, initial->numberOfAtoms()); 
        }
        saddleSearch = new MinModeSaddleSearch(saddle, mode, initial->getPotentialEnergy(), parameters);
    }

    int status = doProcessSearch();

    printEndState(status);
    saveData(status);

    // might have been forced to be equal if the structure passed to the client
    // when determining barrier and the prefactor

    if (min1 != initial) {
        delete initial;
    }
    
    delete saddleSearch;
//    delete initial;
    delete displacement;
    delete saddle;
    delete min1;
    delete min2;

    return returnFiles;
}

int ProcessSearchJob::doProcessSearch(void)
{
    Matter matterTemp(parameters);
    long status;
    int f1;
    f1 = Potential::fcalls;

    
    if (parameters->saddleMethod == "min_mode") {
        status = saddleSearch->run();
    }else{
        DynamicsSaddleSearch dynSearch(initial, parameters);
        int dynSearchStatus = dynSearch.run();

        *saddle = *dynSearch.saddle;
        saddleSearch = dynSearch.saddleSearch;

        if (dynSearchStatus != 0) {
            status = MinModeSaddleSearch::STATUS_BAD_MAX_ITERATIONS;
        }else{
            status = saddleSearch->status;
        }
    }
    fCallsSaddle += Potential::fcalls - f1;

    if (status != MinModeSaddleSearch::STATUS_GOOD) {
        return status;
    }

    // relax from the saddle point

    AtomMatrix posSaddle = saddle->getPositions();
    AtomMatrix displacedPos;

    *min1 = *saddle;
    displacedPos = posSaddle - saddleSearch->getEigenvector() * parameters->processSearchMinimizationOffset;
    min1->setPositions(displacedPos);

    Potential::fcalls = 0;
    log("\nStarting Minimization 1\n");
    bool converged = min1->relax(false, parameters->writeMovies, false, "min1");
    fCallsMin += Potential::fcalls;

    if (!converged) {
        return MinModeSaddleSearch::STATUS_BAD_MINIMA;
    }

    *min2 = *saddle;
    displacedPos = posSaddle + saddleSearch->getEigenvector() * parameters->processSearchMinimizationOffset;
    min2->setPositions(displacedPos);

    Potential::fcalls = 0;
    log("\nStarting Minimization 2\n");
    converged = min2->relax(false, parameters->writeMovies, false, "min2");
    fCallsMin += Potential::fcalls;

    if (!converged) {
        return MinModeSaddleSearch::STATUS_BAD_MINIMA;
    }

    // if min2 corresponds to initial state, swap min1 && min2
    if(!(initial->compare(min1)) && initial->compare(min2)){
        matterTemp = *min1;
        *min1 = *min2;
        *min2 = matterTemp;
    }

    if ((initial->compare(min1)) == false) {
        log("initial != min1\n");
        if((!min1->maxForce() <= parameters->optConvergedForce)  &&
           (!min2->maxForce() <= parameters->optConvergedForce)) {
        }
        return MinModeSaddleSearch::STATUS_BAD_NOT_CONNECTED;
    }

    if (initial->compare(min2)) {
        // both minima are the initial state
        log("both minima are the initial state");
        return MinModeSaddleSearch::STATUS_BAD_NOT_CONNECTED;
    }

    // use the structure passed to the client when determining 
    // the barrier and prefactor for the forward process    
    if (!parameters->processSearchMinimizeFirst) {
        min1 = initial;
    }
    
    // Calculate the barriers
    barriersValues[0] = saddle->getPotentialEnergy()-min1->getPotentialEnergy();
    barriersValues[1] = saddle->getPotentialEnergy()-min2->getPotentialEnergy();

    if((parameters->saddleMaxEnergy < barriersValues[0]) || 
       (parameters->saddleMaxEnergy < barriersValues[1])) {
        return MinModeSaddleSearch::STATUS_BAD_HIGH_BARRIER;
    }

    // calculate the prefactor
    if(!parameters->prefactorDefaultValue)
    {
        f1 = Potential::fcalls;

        int prefStatus;
        double pref1, pref2;
        prefStatus = Prefactor::getPrefactors(parameters, min1, saddle, min2, pref1, pref2);
        if(prefStatus == -1) {
            printf("Prefactor: bad calculation\n");
            return MinModeSaddleSearch::STATUS_FAILED_PREFACTOR;
        }
        fCallsPrefactors += Potential::fcalls - f1;

        /* Check that the prefactors are in the correct range */
        if((pref1 > parameters->prefactorMaxValue) ||
           (pref1 < parameters->prefactorMinValue)){
            cout<<"Bad reactant-to-saddle prefactor: "<<pref1<<endl;
            return MinModeSaddleSearch::STATUS_BAD_PREFACTOR;
        }
        if((pref2 > parameters->prefactorMaxValue) ||
           (pref2 < parameters->prefactorMinValue)){
            cout<<"Bad product-to-saddle prefactor: "<<pref2<<endl;
            return MinModeSaddleSearch::STATUS_BAD_PREFACTOR;
        }
        prefactorsValues[0] = pref1;
        prefactorsValues[1] = pref2;

    }
    else
    {
        // use the default prefactor value specified
        prefactorsValues[0] = parameters->prefactorDefaultValue;
        prefactorsValues[1] = parameters->prefactorDefaultValue;
    }
    return MinModeSaddleSearch::STATUS_GOOD;
}


void ProcessSearchJob::saveData(int status)
{
    FILE *fileResults, *fileReactant, *fileSaddle, *fileProduct, *fileMode;

    std::string resultsFilename("results.dat");
    returnFiles.push_back(resultsFilename);
    fileResults = fopen(resultsFilename.c_str(), "wb");
    ///XXX: min_fcalls isn't quite right it should get them from
    //      the minimizer. But right now the minimizers are in
    //      the SaddleSearch object. They will be taken out eventually.

    fprintf(fileResults, "%d termination_reason\n", status);
    fprintf(fileResults, "%ld random_seed\n", parameters->randomSeed);
    fprintf(fileResults, "%s potential_type\n", parameters->potential.c_str());
    fprintf(fileResults, "%d total_force_calls\n", Potential::fcallsTotal);
    fprintf(fileResults, "%ld force_calls_minimization\n", fCallsMin);
//    fprintf(fileResults, "%ld force_calls_minimization\n", SaddleSearch->forceCallsMinimization + fCallsMin);
    fprintf(fileResults, "%d force_calls_saddle\n", fCallsSaddle);
    fprintf(fileResults, "%f potential_energy_saddle\n", saddle->getPotentialEnergy());
    fprintf(fileResults, "%f potential_energy_reactant\n", min1->getPotentialEnergy());
    fprintf(fileResults, "%f potential_energy_product\n", min2->getPotentialEnergy());
    fprintf(fileResults, "%f barrier_reactant_to_product\n", barriersValues[0]);
    fprintf(fileResults, "%f barrier_product_to_reactant\n", barriersValues[1]);
    if (parameters->saddleMethod == "min_mode") {
        fprintf(fileResults, "%f displacement_saddle_distance\n",
            displacement->perAtomNorm(*saddle));
    }else{
        fprintf(fileResults, "%f displacement_saddle_distance\n", 0.0);
    }
    fprintf(fileResults, "%d force_calls_prefactors\n", fCallsPrefactors);
    fprintf(fileResults, "%.4e prefactor_reactant_to_product\n", prefactorsValues[0]);
    fprintf(fileResults, "%.4e prefactor_product_to_reactant\n", prefactorsValues[1]);
    fclose(fileResults);

    std::string reactantFilename("reactant.con");
    returnFiles.push_back(reactantFilename);
    fileReactant = fopen(reactantFilename.c_str(), "wb");
    min1->matter2con(fileReactant);

    std::string modeFilename("mode_out.dat");
    returnFiles.push_back(modeFilename);
    fileMode = fopen(modeFilename.c_str(), "wb");
    helper_functions::saveMode(fileMode, saddle, saddleSearch->getEigenvector());
    fclose(fileMode);
    fclose(fileReactant);

    std::string saddleFilename("saddle.con");
    returnFiles.push_back(saddleFilename);
    fileSaddle = fopen(saddleFilename.c_str(), "wb");
    saddle->matter2con(fileSaddle);
    fclose(fileSaddle);

    std::string productFilename("product.con");
    returnFiles.push_back(productFilename);
    fileProduct = fopen(productFilename.c_str(), "wb");
    min2->matter2con(fileProduct);
    fclose(fileProduct);

    return;
}

void ProcessSearchJob::printEndState(int status)
{
    fprintf(stdout, "Final state: ");
    if(status == MinModeSaddleSearch::STATUS_GOOD)
        fprintf(stdout, "Successful.\n");

    else if(status == MinModeSaddleSearch::STATUS_BAD_NO_CONVEX)
        fprintf(stdout, "Initial displacement, not able to reach convex region.\n");

    else if(status == MinModeSaddleSearch::STATUS_BAD_HIGH_ENERGY)
        fprintf(stdout, "Saddle search, barrier too high.\n");

    else if(status == MinModeSaddleSearch::STATUS_BAD_MAX_CONCAVE_ITERATIONS)
        fprintf(stdout, "Saddle search, too many iterations in concave region.\n");

    else if(status == MinModeSaddleSearch::STATUS_BAD_MAX_ITERATIONS)
        fprintf(stdout, "Saddle search, too many iterations in saddle point search.\n");

    else if(status == MinModeSaddleSearch::STATUS_BAD_NOT_CONNECTED)
        fprintf(stdout, "Minima, saddle is not connected to initial state.\n");

    else if(status == MinModeSaddleSearch::STATUS_BAD_PREFACTOR)
            fprintf(stdout, "Prefactors, not within window\n");

    else if(status == MinModeSaddleSearch::STATUS_FAILED_PREFACTOR)
            fprintf(stdout, "Prefactors, hessian calculation failed\n");

    else if(status == MinModeSaddleSearch::STATUS_BAD_HIGH_BARRIER)
        fprintf(stdout, "Energy barriers, not within window\n");

    else if (status == MinModeSaddleSearch::STATUS_BAD_MINIMA)
        fprintf(stdout, "Minima, from saddle did not converge\n");

    else if(status == MinModeSaddleSearch::STATUS_NONNEGATIVE_ABORT)
        log("[SaddleSearch] Nonnegative initial mode, aborting.\n");

    else
        fprintf(stdout, "Unknown status: %i!\n", status);

    return;
}

