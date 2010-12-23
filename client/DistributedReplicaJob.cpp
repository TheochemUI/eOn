//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include <cstdlib>

#include "Matter.h"
#include "Dynamics.h"
#include "DistributedReplicaJob.h"
#include "ConjugateGradients.h"

#ifdef BOINC
    #include <boinc/boinc_api.h>
    #include <boinc/diagnostics.h>    
    #include <boinc/filesys.h>        
#ifdef WIN32
    #include <boinc/boinc_win.h>
    #include <boinc/win_util.h>
#endif
#else
    #include "false_boinc.h"
#endif

DistributedReplicaJob::DistributedReplicaJob(Parameters *params)
{
    double target_temp;
    parameters = params;
    min_fcalls = 0;
    bl_fcalls = 0;
    sp_fcalls = 0;
    rf_fcalls = 0;
    temperature = parameters->temperature;
    target_temp = parameters->drTargetTemperature;
    save_refine = false;
    if (target_temp == temperature){
        save_refine = true;
    }
//    printf("%lf, %lf, %d",temperature,target_temp,save_refine);
}

DistributedReplicaJob::~DistributedReplicaJob()
{ 

}

void DistributedReplicaJob::run(int bundleNumber)
{
    char buff[STRING_SIZE];
    string reactant_passed("reactant_passed");

    if (bundleNumber < 0) {
        reactant_passed += ".con";
    }else{
        snprintf(buff, STRING_SIZE, "_%i.con", bundleNumber);
        reactant_passed += buff;
    }

    reactant = new Matter(parameters);
    min1 = new Matter(parameters);
    min2 = new Matter(parameters);
    final = new Matter(parameters);

    reactant->convel2matter(reactant_passed);

  /*
    Matrix<double, Eigen::Dynamic, 3> vel;
    vel = reactant->getVelocities();
    printf("%lf %lf %lf\n",vel(0,0),vel(0,1),vel(0,2));
 */

    *min1 = *reactant;
    *min2 = *reactant;
    *final = *reactant;

    ConjugateGradients cgMin1(min1, parameters);
    cgMin1.setOutput(0);
    printf("\nMinimizing initial reactant\n");
    cgMin1.fullRelax();
    min_fcalls += min1->getForceCalls();

    printf("Now running Distributed Replica Simulation\n\n"); 

    balanceStep();
    samplingStep();
    saveData(bundleNumber);

    delete reactant;
    delete min1;
    delete min2;
    delete final;
}

void DistributedReplicaJob::balanceStep(){
    long n, bSteps,nloop;
    bool bl_new,stop;
    Matter *initial;
    
    n = 0;
    nloop = 0;
    bl_new = false;
    stop = false;
    bSteps = parameters->drBalanceSteps;
    initial = new Matter(parameters);
    *initial = *reactant;

    Dynamics balanceDynamics(reactant,parameters);
    printf("Balancing for %ld steps :",bSteps);

    balanceDynamics.velRescaling(parameters->temperature);
    while(!stop){
        while(n < bSteps){
           balanceDynamics.oneStep(parameters->temperature);
           n++;
           bl_fcalls++;
        }
        bl_new = checkState(reactant);
        if(bl_new){
            *reactant = *initial;
            stop = false;
            nloop ++;
        }else{
            printf(" Succeed!\n");
            stop = true;
        }
        if(nloop >= 5){ 
            *reactant = *initial;
            printf("Warning - The system is in a meta-stable state!\n");
            stop = true;
        }
    }
    return;
    delete initial;
}

void DistributedReplicaJob::samplingStep(){
    long n, sSteps;
    n = 0;
    sSteps = parameters->drBalanceSteps;

    Dynamics samplingDynamics(reactant,parameters);
    printf("Sampling for %ld steps\n",sSteps);

    while(n < sSteps){
        samplingDynamics.oneStep(parameters->temperature);
        n++;
        sp_fcalls++;
    }

    *final = *reactant;
    *min2 = *final;
    ConjugateGradients cgMin2(min2, parameters);     
    cgMin2.setOutput(0);     
    printf("\nMinimizing final product\n");     
    cgMin2.fullRelax();     
    min_fcalls += min2->getForceCalls();

    return;
}

void DistributedReplicaJob::saveData(int bundleNumber){
 
    FILE *fileResults, *fileReactant, *fileProduct;

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
    long total_fcalls = min_fcalls + bl_fcalls + sp_fcalls + rf_fcalls;

    fprintf(fileResults, "%ld random_seed\n", parameters->randomSeed);
    fprintf(fileResults, "%lf potential_energy_reactant\n", min1->getPotentialEnergy());
    fprintf(fileResults, "%lf potential_energy_product\n", final->getPotentialEnergy());
    fprintf(fileResults, "%ld potential_type\n", parameters->potential);
    fprintf(fileResults, "%ld force_calls_balance\n", bl_fcalls);
    fprintf(fileResults, "%ld force_calls_sampling\n", sp_fcalls);
    fprintf(fileResults, "%ld force_calls_minimization\n", min_fcalls);
    fprintf(fileResults, "%ld force_calls_total\n", total_fcalls);

    fprintf(fileResults, "%lf moved_distance\n",min2->distanceTo(*min1));
    fclose(fileResults);

    if (bundleNumber != -1) {
        snprintf(filename, STRING_SIZE, "reactant_%i.con", bundleNumber);
    }else{
         strncpy(filename, "reactant.con", STRING_SIZE);
    }
    fileReactant = fopen(filename, "wb");
    min1->matter2con(fileReactant);
    fclose(fileReactant);

    if (bundleNumber != -1) {
        snprintf(filename, STRING_SIZE, "product_%i.convel", bundleNumber);
    }else{
        strncpy(filename, "product.convel", STRING_SIZE);
    }

    fileProduct = fopen(filename, "wb");
    final->matter2convel(fileProduct);
    fclose(fileProduct);
 
    return;

}


bool DistributedReplicaJob::checkState(Matter *matter)
{
    Matter tmp(parameters);
    tmp = *matter;
    ConjugateGradients cgMin(&tmp, parameters);
    cgMin.fullRelax();
    min_fcalls += tmp.getForceCalls();

    if (tmp == *min1) {
        return false;
    }
    return true;
}

