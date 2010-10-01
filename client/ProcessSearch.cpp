#include "ProcessSearch.h"
#include "Constants.h"
#include "false_boinc.h"

#include <stdio.h>
#include <string>

using namespace std;

ProcessSearch::ProcessSearch (Parameters *params)
{
    parameters = params;
}

ProcessSearch::~ProcessSearch()
{

}

void ProcessSearch::run(int bundleNumber)
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

    if (parameters->getRefineSP()) {
        saddle->con2matter(displacement_passed);
    }

    barriersValues[0] = barriersValues[1] = 0;
    prefactorsValues[0] = prefactorsValues[1] = 0;

    if (parameters->getRefineSP()) {
        *min1 = *min2 = *initial;
    }else{
        *saddle = *min1 = *min2 = *initial;
    }

    saddlePoint = new SaddlePoint();
    saddlePoint->initialize(initial, saddle, parameters);
    if (parameters->getRefineSP()) {
        saddlePoint->loadMode(mode_passed);
    }

    prefactors = new Prefactors();
    prefactors->initialize(saddle, min1, min2, parameters);

    int status = doProcessSearch();

    printEndState(status);
    saveData(status, bundleNumber);

    delete prefactors;
    delete saddlePoint;
    delete initial;
    delete displacement;
    delete saddle; 
    delete min1; 
    delete min2; 
}

int ProcessSearch::doProcessSearch(void)
{
    Matter matterTemp(parameters);
    long status;
    status = saddlePoint->locate(min1, min2);

    if (status != statusInit) {
        return status;
    }

    // If min2 corresponds to initial state swap min1 && min2
    if(!(*initial==*min1) && ((*initial==*min2))){
        matterTemp = *min1;
        *min1 = *min2;
        *min2 = matterTemp;
    }

    if ((*initial==*min1) == false) {
        printf("initial != min1\n");
        if((!min1->isItConverged(parameters->getConverged_Relax()))  &&
           (!min2->isItConverged(parameters->getConverged_Relax()))) {
            return statusBadMinima;
        }
        return statusBadNotConnected;
    }

    //if((!min1->isItConverged(parameters->getConverged_Relax())))
    //{
    //    printf("min1 is not converged!!!!!! %lf\n", parameters->getConverged_Relax());
    //}
    //if((!min2->isItConverged(parameters->getConverged_Relax())))
    //{
    //    printf("min2 is not converged!!!!!! %lf\n", parameters->getConverged_Relax());
    //}
    
    
    
    if (*initial==*min2) {
        /* both minima are the initial state */
        printf("both minima are the initial state");
        return statusBadNotConnected;
    }


    /* Calculate the barriers */
    barriersValues[0] = saddle->potentialEnergy()-min1->potentialEnergy();
    barriersValues[1] = saddle->potentialEnergy()-min2->potentialEnergy();

    if((parameters->getMaxEnergy_SP() < barriersValues[0]) || 
       (parameters->getMaxEnergy_SP() < barriersValues[1])) {
        return statusBadHighBarrier;
    }

    /* Perform the dynamical matrix caluclation */
    prefactors->compute(prefactorsValues);

    /* Check that the prefactors are in the correct range */
    if((prefactorsValues[0]>parameters->getPrefactorMax()) ||
       (prefactorsValues[0]<parameters->getPrefactorMin())){
        return statusBadPrefactor;
    }
    if((prefactorsValues[1]>parameters->getPrefactorMax()) ||
       (prefactorsValues[1]<parameters->getPrefactorMin())){
        return statusBadPrefactor;
    }

    return statusGood;
}

void ProcessSearch::saveData(int status, int bundleNumber){
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
    long min_fcalls = min1->getForceCalls()+min2->getForceCalls();
    long saddle_fcalls = saddlePoint->forceCallsSaddlePointConcave_ + 
                         saddlePoint->forceCallsSaddlePointConvex_;
    long total_fcalls = min_fcalls + saddle_fcalls+prefactors->totalForceCalls;

    fprintf(fileResults, "%d termination_reason\n", status);
    fprintf(fileResults, "%ld random_seed\n", parameters->getRandomSeed());
    fprintf(fileResults, "%ld potential_tag\n", parameters->getPotentialTag());
    fprintf(fileResults, "%ld total_force_calls\n", total_fcalls);
    fprintf(fileResults, "%ld force_calls_minimization\n", min_fcalls);
    fprintf(fileResults, "%ld force_calls_saddle\n", saddle_fcalls);
    fprintf(fileResults, "%f potential_energy_saddle\n", saddle->potentialEnergy());
    fprintf(fileResults, "%f potential_energy_reactant\n", min1->potentialEnergy());
    fprintf(fileResults, "%f potential_energy_product\n", min2->potentialEnergy());
    fprintf(fileResults, "%f barrier_reactant_to_product\n", barriersValues[0]);
    fprintf(fileResults, "%f barrier_product_to_reactant\n", barriersValues[1]);
    fprintf(fileResults, "%f displacement_saddle_distance\n",
            displacement->per_atom_norm(*saddle));
    fprintf(fileResults, "%ld force_calls_prefactors\n", prefactors->totalForceCalls);
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

void ProcessSearch::printEndState(int status) {
    fprintf(stdout, "Final state: ");
    if(status == statusGood)
        fprintf(stdout, "Successful.\n");

    else if(status == statusBadNoConvex)
        fprintf(stdout, "Initial displacement, not able to reach convex region.\n");
   
    else if(status == statusBadHighEnergy)
        fprintf(stdout, "Saddle search, barrier too high.\n");
        
    else if(status == statusBadMaxConcaveIterations) 
        fprintf(stdout, "Saddle search, too many iterations in concave region.\n");

    else if(status == statusBadMaxIterations)
        fprintf(stdout, "Saddle search, too many iterations in saddle point search.\n");

    else if(status == statusBadNotConnected)
        fprintf(stdout, "Minima, saddle is not connected to initial state.\n");

    else if(status == statusBadPrefactor)
            fprintf(stdout, "Prefactors, not within window as defined in Constants\n");

    else if(status == statusBadHighBarrier)
        fprintf(stdout, "Energy barriers, not within window as defined in Constants\n");
    else if (status == statusBadMinima)
        fprintf(stdout, "Minima were not able to be matched to the initial reactant\n"
                "because they did not converge\n");
    else
        fprintf(stdout, "Unknown status: %i!\n", status);

    return;
}

