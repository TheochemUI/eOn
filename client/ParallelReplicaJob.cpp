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
#include "BondBoost.h"
#include "ParallelReplicaJob.h"
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

ParallelReplicaJob::ParallelReplicaJob(Parameters *params)
{
    parameters = params;
    nsteps = 0;
    nsteps_refined = 0;
    SPtime = 0.0;
    RLtime = 0.0;
    check_steps = parameters->mdCheckFreq;
    relax_steps = parameters->mdRelaxSteps;
    SPtimebuff = new double[check_steps];
    newstate = false;
    meta = false;
    min_fcalls = 0;
    md_fcalls = 0;
    rf_fcalls = 0;
    dh_fcalls = 0;
}

ParallelReplicaJob::~ParallelReplicaJob()
{ 
    delete[] SPtimebuff;
}

void ParallelReplicaJob::run(int bundleNumber)
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
    fin1 = new Matter(parameters);
    fin2 = new Matter(parameters);
    saddle = new Matter(parameters);
    final = new Matter(parameters);
    
    reactant->con2matter(reactant_passed);
    *min1 = *reactant;   
    *saddle = *reactant;
    *fin1 = *reactant;
    *fin2 = *reactant;

    ConjugateGradients cgMin1(min1, parameters);
    cgMin1.setOutput(0);
    printf("\nMinimizing initial reactant\n");
    cgMin1.fullRelax();
    min_fcalls += min1->getForceCalls();
    *final = *min1;    

    printf("Now running Parallel Replica Dynamics\n\n");

    dynamics();
    saveData(newstate,bundleNumber);
    
    if(newstate){
     //   printf("New state has been found\n");
        printf("Transition Time: %.2e s\n", SPtime*1.018e-14);
    }else{
       printf("New state has not been found in this %ld dynamics steps (%.2f fs)\n",
            parameters->mdSteps,10.18*parameters->mdSteps*parameters->mdTimeStep);
    }

    delete min1;
    delete fin1;
    delete fin2;
    delete reactant;
    delete saddle;
    delete final;
}

void ParallelReplicaJob::dynamics()
{
    bool   status = false, remember = true, stoped = false;
    long   nFreeCoord = reactant->numberOfFreeAtoms()*3;
    long   RecordAccuracy = parameters->mdRecordAccuracy, mdbufflength;
    long   ncheck = 0, nexam = 0, nrecord = 0, steps_tmp = 0, final_refined;
    Matrix<double, Eigen::Dynamic, 3> velocities;
    double kinE = 0.0, kb = 1.0/11604.5;
    double kinT = 0.0, sumT = 0.0, sumT2 = 0.0, avgT, varT;

    mdbufflength = int(check_steps/RecordAccuracy)+1;
    Matter *mdbuff[mdbufflength];
    for(long i =0; i < mdbufflength;i++){
        mdbuff[i] = new Matter(parameters);
    }
    printf("RecordAccuracy = %ld; mdbufflength = %ld\n",RecordAccuracy,mdbufflength);
    
    Dynamics PRdynamics(reactant,parameters);
    BondBoost Bbm(reactant,parameters);

    if(parameters->biasPotential == Hyperdynamics::BOND_BOOST){
        Bbm.initial();
    }

    PRdynamics.initialVel(parameters->temperature);
    dephase();

    printf("\nStarting MD run\nTemperature: %.2f Kelvin\nTotal Time: %.2f fs\nTime Step: %.2f fs\n\n",
    parameters->temperature, 10.18*parameters->mdSteps*parameters->mdTimeStep, 10.18*parameters->mdTimeStep);

    long tenthSteps = parameters->mdSteps/10;
    //This prevents and edge case division by zero if mdSteps is < 10
    if (tenthSteps == 0) {
        tenthSteps = parameters->mdSteps;
    }
   
    while(!stoped){

        if(parameters->biasPotential == Hyperdynamics::BOND_BOOST && !newstate){
            SPtime += Bbm.boost();
        }

        else{ SPtime += parameters->mdTimeStep;}
                
        kinE = reactant->getKineticEnergy();
        kinT = (2*kinE/nFreeCoord/kb); 
        sumT += kinT;
        sumT2 += kinT*kinT;

        PRdynamics.oneStep(parameters->temperature);

        ncheck++;
        nsteps++;

        if(parameters->mdRefine && remember && !newstate ){
            if(ncheck % RecordAccuracy == 0){
                nrecord ++;
                *mdbuff[nrecord-1] = *reactant;
                SPtimebuff[nrecord-1] = SPtime;
            }
        }
        
       // printf("MDsteps %ld Ekin= %lf Tkin= %lf \n",nsteps,kinE,kinT); 

#ifndef NDEBUG
        if (ncheck == check_steps && !newstate){
           reactant->matter2xyz("movie", true);
        }
#endif
        if (ncheck == check_steps && !newstate){
            ncheck = 0; // reinitialize the ncheck
            nrecord = 0;
            status = PRdynamics.checkState(reactant,min1);
            if(status == true){
                remember = false;
            }
        }

        if (status && !newstate){
            nexam++;
            if (nexam >= relax_steps){
                nexam = 0;
                ncheck = 0; 	
                nrecord = 0;
                newstate = PRdynamics.checkState(reactant,min1);
                status = false;
                if(newstate == false){
                    remember = true;
                }else{
                    printf("Found New State !\n");
                    *final = * reactant;
                    steps_tmp = nsteps;
                    if(parameters->mdAutoStop){
                        stoped = true;
                    }
                    remember = false;
                }
            }
        }

        if (nsteps >= parameters->mdSteps){
           stoped = true;
        }
 
        //BOINC Progress
        if (nsteps % 500 == 0) {
            // Since we only have a bundle size of 1 we can play with boinc_fraction_done
            // directly. When we have done parameters->mdSteps number of steps we aren't 
            // quite done so I increase the max steps by 5%.
            boinc_fraction_done((double)nsteps/(double)(parameters->mdSteps+0.05*parameters->mdSteps));
        }

        //stdout Progress
        if (nsteps % tenthSteps == 0 || nsteps == parameters->mdSteps) {
            double maxAtomDistance = reactant->perAtomNorm(*min1);
            printf("progress: %3.0f%%, max displacement: %6.3lf, step %7ld/%ld\n",
                   (double)100.0*nsteps/parameters->mdSteps, maxAtomDistance,
                   nsteps, parameters->mdSteps);
        }
    }
    avgT=sumT/nsteps;
    varT=sumT2/nsteps-avgT*avgT;
    printf("\nTemperature : Average = %lf ; Variance = %lf ; Factor = %lf\n\n", avgT,varT,varT/avgT/avgT*nFreeCoord/2);

    if (isfinite(avgT)==0) {
        printf("Infinite average temperature, something went wrong!\n");
        newstate = false;
    }

    final_refined = steps_tmp;

    //printf("mdbufflength = %ld; nrecord = %ld\n",mdbufflength, nrecord);
    if(parameters->mdRefine && newstate){
        nsteps_refined = PRdynamics.refine(mdbuff,mdbufflength,min1);
        //printf("nsteps_refined=%ld\n",nsteps_refined); 
        final_refined = steps_tmp-check_steps-relax_steps+nsteps_refined*RecordAccuracy;
        //long totsteps = nsteps-check_steps-relax_steps+nsteps_refined+1; 
        //printf("final_step = %ld\n",final_refined);
        *reactant = *mdbuff[nsteps_refined];
        *saddle = *reactant;
        SPtime = SPtimebuff[nsteps_refined];
    
        printf("Found transition at step %ld, now running another %ld steps to allocate the product state\n",final_refined, relax_steps);
        
        long relaxbufflength = int(relax_steps/RecordAccuracy)+1;
        if(nsteps_refined < mdbufflength - relaxbufflength ){
            *final = *mdbuff[nsteps_refined+relaxbufflength];
        }
    
        *fin1 = *saddle;
        ConjugateGradients SaddleMin(fin1, parameters);
        SaddleMin.fullRelax();
        min_fcalls += fin1->getForceCalls();


        *fin2 = *final;
        ConjugateGradients RelaxMin(fin2, parameters);
        RelaxMin.fullRelax();
        min_fcalls += fin2->getForceCalls();

        if(*fin2 == *fin1){
           *final = *fin1;
           printf("Transition followed by a stable state !\n");
        }else{
           *final = *fin2;
           printf("Transition followed by a meta-stable state; using fin2 as product.con\n");
           SPtime = SPtime + relax_steps;
           meta = true;
        }
    }

    min_fcalls += PRdynamics.getMinfcalls();
    rf_fcalls += PRdynamics.getRefinefcalls();
    md_fcalls += PRdynamics.getMDfcalls();

    return;

    for(long i =0; i < mdbufflength;i++){
        delete[] mdbuff[i];
    }
}

void ParallelReplicaJob::saveData(int status,int bundleNumber)
{
    FILE *fileResults, *fileReactant, *fileProduct, *fileSaddle, *fileMeta;

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
    long total_fcalls = min_fcalls + md_fcalls + dh_fcalls + rf_fcalls;

    fprintf(fileResults, "%d termination_reason\n", status);
    fprintf(fileResults, "%e transition_time_s\n", SPtime*1.018e-14);
    fprintf(fileResults, "%ld random_seed\n", parameters->randomSeed);
    fprintf(fileResults, "%lf potential_energy_reactant\n", min1->getPotentialEnergy());
    fprintf(fileResults, "%lf potential_energy_product\n", final->getPotentialEnergy());
    fprintf(fileResults, "%ld potential_type\n", parameters->potential);
    fprintf(fileResults, "%ld total_force_calls\n", total_fcalls);
    fprintf(fileResults, "%ld force_calls_dephase\n", dh_fcalls);
    fprintf(fileResults, "%ld force_calls_dynamics\n", md_fcalls);
    fprintf(fileResults, "%ld force_calls_minimization\n", min_fcalls);
    fprintf(fileResults, "%ld force_calls_refine\n", rf_fcalls);

    fprintf(fileResults, "%lf moved_distance\n",final->distanceTo(*min1));
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
        snprintf(filename, STRING_SIZE, "product_%i.con", bundleNumber);
    }else{
        strncpy(filename, "product.con", STRING_SIZE);
    }

    fileProduct = fopen(filename, "wb");
    final->matter2con(fileProduct);
    fclose(fileProduct);
  
    if (bundleNumber != -1) {
        snprintf(filename, STRING_SIZE, "saddle_%i.con", bundleNumber);
    }else{
        strncpy(filename, "saddle.con", STRING_SIZE);
    }

    fileSaddle = fopen(filename, "wb");
    saddle->matter2con(fileSaddle);
    fclose(fileSaddle);

    if(meta){
       if (bundleNumber != -1) {
           snprintf(filename, STRING_SIZE, "meta_%i.con", bundleNumber);
       }else{
           strncpy(filename, "meta.con", STRING_SIZE);
       }

       fileMeta = fopen(filename, "wb");
       fin1->matter2con(fileMeta);
       fclose(fileMeta);
    }
    return;
}

void ParallelReplicaJob::dephase()
{
    long  dephaseSteps = parameters->mdDephaseSteps,i = 0;
    long  n,nloop, new_n, nbuff, dh_refined;
    bool  state = false;
    Matrix<double, Eigen::Dynamic, 3> velocity;      

    Dynamics dephaseDynamics(reactant,parameters);
    printf("Dephasing for %ld steps\n",dephaseSteps);
    
    n  = 0;
    new_n = 0;
    nloop= 0;
    while(n < dephaseSteps){

        nbuff = dephaseSteps-n;
        nloop ++;
        Matter *DHbuff[nbuff];
        for(i=0; i<nbuff; i++){
           DHbuff[i] = new Matter(parameters);
        }

        for(i=0; i<nbuff; i++){
           dephaseDynamics.oneStep(parameters->temperature);
           *DHbuff[i] = *reactant;
        }
        
        state = dephaseDynamics.checkState(reactant,min1);
        if(state == true){
            dh_refined = dephaseDynamics.refine(DHbuff,nbuff,min1);
            printf("nloop = %ld; refined step = %ld\n",nloop,dh_refined);
            new_n = dh_refined-parameters->mdRefineAccuracy;
            new_n = (new_n > 0)?new_n:0;
            printf("Dephasing Warning: in a new state,now inverse the momentum and restart from step %ld\n",n+new_n);
            *reactant = *DHbuff[new_n];
            velocity = reactant->getVelocities();
            velocity = velocity*(-1);
            reactant->setVelocities(velocity);
            n = n + new_n;
        }else{
            n = n + nbuff;
            printf("Successful Dephasing for %ld steps \n",n);            
        }
   
       for(i=0; i<nbuff;i++){
           delete DHbuff[i];
       }

       if(parameters->mdDephaseLoopStop){
           if(nloop > parameters->mdDephaseLoopMax){
              printf("Reach Dephase Loop Maximum, Stop Dephasing! Dephsed for %ld steps\n ",n);
              break;    
           }
       }
   }
   min_fcalls += dephaseDynamics.getMinfcalls(); 
   rf_fcalls += dephaseDynamics.getRefinefcalls();
   dh_fcalls += dephaseDynamics.getMDfcalls();
}

