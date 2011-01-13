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
#include "Constants.h"
#include "ConjugateGradients.h"
#include "false_boinc.h"
#include "Potentials.h"

#include <stdio.h>
#include <string>

using namespace std;

ProcessSearchJob::ProcessSearchJob (Parameters *params)
{
    parameters = params;
    fCallsSaddle = fCallsPrefactors = fCallsMin = 0;
}

ProcessSearchJob::~ProcessSearchJob()
{}

void ProcessSearchJob::run(int bundleNumber)
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
    min1 = new Matter(parameters);
    min2 = new Matter(parameters);

    initial->con2matter(reactant_passed);

    if (parameters->processSearchMinimizeFirst) {
        printf("Minimizing initial structure\n");
        int fi = Potential::fcalls;
        ConjugateGradients cgMin(initial, parameters);
        cgMin.fullRelax();
        fCallsMin += Potential::fcalls - fi;
    }

    if (!parameters->saddleDisplaceType) {
        saddle->con2matter(displacement_passed);
    }
    
    barriersValues[0] = barriersValues[1] = 0;
    prefactorsValues[0] = prefactorsValues[1] = 0;

    if (!parameters->saddleDisplaceType) {
        *min1 = *min2 = *initial;
    }else{
        *saddle = *min1 = *min2 = *initial;
    }

    saddlePoint = new SaddlePoint();
    saddlePoint->initialize(initial, saddle, parameters);
    if (!parameters->saddleDisplaceType) {
        saddlePoint->loadMode(mode_passed);
    }
    else{
        saddlePoint->displaceAndSetMode(saddle);
    }
    

    hessian = new Hessian(min1, saddle, min2, parameters);

    int status = doProcessSearch();

    printEndState(status);
    saveData(status, bundleNumber);

    delete hessian;
    delete saddlePoint;
    delete initial;
    delete displacement;
    delete saddle; 
    delete min1; 
    delete min2; 
}

int ProcessSearchJob::doProcessSearch(void)
{
    Matter matterTemp(parameters);
    long status;
    int f1;
    f1 = Potential::fcalls;
    status = saddlePoint->locate();//(min1, min2);
    fCallsSaddle += Potential::fcalls - f1;

    if (status != SaddlePoint::STATUS_INIT) {
        return status;
    }
	
	// relax from the saddle point located
	
	Matrix<double, Eigen::Dynamic, 3> posSaddle = saddlePoint->getSaddlePositions();
	Matrix<double, Eigen::Dynamic, 3> displacedPos;

	*min1 = *saddle;
	//XXX: the distance displaced from the saddle should be a parameter
	displacedPos = posSaddle - saddlePoint->getEigenMode() * 0.2;
	min1->setPositions(displacedPos);
	ConjugateGradients cgMin1(min1, parameters);
	cgMin1.fullRelax();
	fCallsMin += cgMin1.totalForceCalls;
		
	*min2 = *saddle;
	displacedPos = posSaddle + saddlePoint->getEigenMode() * 0.2;
	min2->setPositions(displacedPos);
	ConjugateGradients cgMin2(min2, parameters);  
	cgMin2.fullRelax();
	fCallsMin += cgMin2.totalForceCalls;
	
    // If min2 corresponds to initial state swap min1 && min2
    if(!(*initial==*min1) && ((*initial==*min2))){
        matterTemp = *min1;
        *min1 = *min2;
        *min2 = matterTemp;
    }

    if ((*initial==*min1) == false) {
        printf("initial != min1\n");
        if((!min1->isItConverged(parameters->optConvergedForce))  &&
           (!min2->isItConverged(parameters->optConvergedForce))) {
            // the isItConverged in Matter doesn't work so this needs to be fixed.
            // return statusBadMinima;
        }
        return SaddlePoint::STATUS_BAD_NOT_CONNECTED;
    }

    //if((!min1->isItConverged(parameters->convergedRelax)))
    //{
    //    printf("min1 is not converged!!!!!! %lf\n", parameters->convergedRelax);
    //}
    //if((!min2->isItConverged(parameters->convergedRelax)))
    //{
    //    printf("min2 is not converged!!!!!! %lf\n", parameters->convergedRelax);
    //}

    if (*initial==*min2) {
        /* both minima are the initial state */
        printf("both minima are the initial state");
        return SaddlePoint::STATUS_BAD_NOT_CONNECTED;
    }

    /* Calculate the barriers */
    barriersValues[0] = saddle->getPotentialEnergy()-min1->getPotentialEnergy();
    barriersValues[1] = saddle->getPotentialEnergy()-min2->getPotentialEnergy();

    if((parameters->saddleMaxEnergy < barriersValues[0]) || 
       (parameters->saddleMaxEnergy < barriersValues[1])) {
        return SaddlePoint::STATUS_BAD_HIGH_BARRIER;
    }

    if(!parameters->processSearchDefaultPrefactor)
    {
        /* Perform the dynamical matrix caluclation */
        VectorXd reactModes, saddleModes, prodModes;
        f1 = Potential::fcalls;
        reactModes = hessian->getModes(Hessian::REACTANT);
        if(reactModes.size() == 0)
        {
            return SaddlePoint::STATUS_BAD_PREFACTOR;
        }
        saddleModes = hessian->getModes(Hessian::SADDLE);
        if(saddleModes.size() == 0)
        {
            return SaddlePoint::STATUS_BAD_PREFACTOR;
        }
        prodModes = hessian->getModes(Hessian::PRODUCT);
        if(prodModes.size() == 0)
        {
            return SaddlePoint::STATUS_BAD_PREFACTOR;
        }
        fCallsPrefactors += Potential::fcalls - f1; 

        prefactorsValues[0] = 1;
        prefactorsValues[1] = 1;
    
        // This may have large numerical error
        // products are calculated this way in order to avoid overflow
        for(int i=0; i<saddleModes.size(); i++)
        {
            prefactorsValues[0] *= reactModes[i];
            prefactorsValues[1] *= prodModes[i];
            if(saddleModes[i]>0)
            {
                prefactorsValues[0] /= saddleModes[i];
                prefactorsValues[1] /= saddleModes[i];
            }

        }
        prefactorsValues[0] = sqrt(prefactorsValues[0])/(2*M_PI*10.18e-15);
        prefactorsValues[1] = sqrt(prefactorsValues[1])/(2*M_PI*10.18e-15);

        /* Check that the prefactors are in the correct range */
        if((prefactorsValues[0]>parameters->hessianPrefactorMax) ||
           (prefactorsValues[0]<parameters->hessianPrefactorMin)){
            return SaddlePoint::STATUS_BAD_PREFACTOR;
        }
        if((prefactorsValues[1]>parameters->hessianPrefactorMax) ||
           (prefactorsValues[1]<parameters->hessianPrefactorMin)){
            return SaddlePoint::STATUS_BAD_PREFACTOR;
        }
    }else{ // use the default prefactor value specified
        prefactorsValues[0]=parameters->processSearchDefaultPrefactor;
        prefactorsValues[1]=parameters->processSearchDefaultPrefactor;
    }
    return SaddlePoint::STATUS_GOOD;
}

void ProcessSearchJob::saveData(int status, int bundleNumber){
    FILE *fileResults, *fileReactant, *fileSaddle, *fileProduct, *fileMode;

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
    fprintf(fileResults, "%ld potential_type\n", parameters->potential);
    fprintf(fileResults, "%d total_force_calls\n", Potential::fcalls);
    fprintf(fileResults, "%ld force_calls_minimization\n", fCallsMin);
//    fprintf(fileResults, "%ld force_calls_minimization\n", saddlePoint->forceCallsMinimization + fCallsMin);
    fprintf(fileResults, "%d force_calls_saddle\n", fCallsSaddle);
    fprintf(fileResults, "%f potential_energy_saddle\n", saddle->getPotentialEnergy());
    fprintf(fileResults, "%f potential_energy_reactant\n", min1->getPotentialEnergy());
    fprintf(fileResults, "%f potential_energy_product\n", min2->getPotentialEnergy());
    fprintf(fileResults, "%f barrier_reactant_to_product\n", barriersValues[0]);
    fprintf(fileResults, "%f barrier_product_to_reactant\n", barriersValues[1]);
    fprintf(fileResults, "%f displacement_saddle_distance\n",
            displacement->perAtomNorm(*saddle));
    fprintf(fileResults, "%d force_calls_prefactors\n", fCallsPrefactors);
    fprintf(fileResults, "%.4e prefactor_reactant_to_product\n", prefactorsValues[0]);
    fprintf(fileResults, "%.4e prefactor_product_to_reactant\n", prefactorsValues[1]);
	fclose(fileResults);

    if (bundleNumber != -1) {
        snprintf(filename, STRING_SIZE, "reactant_%i.con", bundleNumber);
    }else{
        strncpy(filename, "reactant.con", STRING_SIZE);
    }
	fileReactant = fopen(filename, "wb");
    min1->matter2con(fileReactant);

    if (bundleNumber != -1) {
        snprintf(filename, STRING_SIZE, "mode_%i.dat", bundleNumber);
    }else{
        strncpy(filename, "mode.dat", STRING_SIZE);
    }
    fileMode = fopen(filename, "wb");
    saddlePoint->saveMode(fileMode);
    fclose(fileMode);
	fclose(fileReactant);

    if (bundleNumber != -1) {
        snprintf(filename, STRING_SIZE, "saddle_%i.con", bundleNumber);
    }else{
        strncpy(filename, "saddle.con", STRING_SIZE);
    }
	fileSaddle = fopen(filename, "wb");
    saddle->matter2con(fileSaddle);
	fclose(fileSaddle);

    if (bundleNumber != -1) {
        snprintf(filename, STRING_SIZE, "product_%i.con", bundleNumber);
    }else{
        strncpy(filename, "product.con", STRING_SIZE);
    }
	fileProduct = fopen(filename, "wb");
	min2->matter2con(fileProduct);
    fclose(fileProduct);

    return;
}

void ProcessSearchJob::printEndState(int status) {
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

    else if(status == SaddlePoint::STATUS_BAD_NOT_CONNECTED)
        fprintf(stdout, "Minima, saddle is not connected to initial state.\n");

    else if(status == SaddlePoint::STATUS_BAD_PREFACTOR)
            fprintf(stdout, "Prefactors, not within window\n");

    else if(status == SaddlePoint::STATUS_BAD_HIGH_BARRIER)
        fprintf(stdout, "Energy barriers, not within window\n");

    else if (status == SaddlePoint::STATUS_BAD_MINIMA)
        fprintf(stdout, "Minima, from saddle did not converge\n");
    else
        fprintf(stdout, "Unknown status: %i!\n", status);

    return;
}

