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
    min_fcalls = 0;
    md_fcalls = 0;
    rf_fcalls = 0;
    dh_fcalls = 0;
}

ParallelReplicaJob::~ParallelReplicaJob(){ }

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
    min2 = new Matter(parameters);
    transition = new Matter(parameters);

    reactant->con2matter(reactant_passed);
    *min1 = *reactant;   
    *transition = *reactant;
    ConjugateGradients cgMin1(min1, parameters);
    cgMin1.setOutput(0);
    printf("\nMinimizing initial reactant\n");
    cgMin1.fullRelax();
    min_fcalls += min1->getForceCalls();
        
    printf("Now running Parallel Replica Dynamics\n\n");

    dynamics();

    *min2 = *transition;   
    ConjugateGradients cgMin2(min2, parameters);
    cgMin2.fullRelax();
    min_fcalls += min2->getForceCalls();

    saveData(newstate,bundleNumber);
    
    if(newstate){
        printf("New state has been found\n");
        printf("Transition Time: %.2e\n", SPtime*1e-15);
    }else{
       printf("New state has not been found in this %ld dynamics steps (%.2f fs)\n",
            parameters->mdSteps,10.18*parameters->mdSteps*parameters->mdTimeStep);
    }

    delete min1;
    delete min2;
    delete reactant;
    delete transition;
    delete[] SPtimebuff;
}

void ParallelReplicaJob::dynamics()
{
    bool   status = false, remember = true, stoped = false;
    bool   boost = parameters->bondBoost;
    long   nFreeCoord = reactant->numberOfFreeAtoms()*3;
    long   ncheck = 0, nexam = 0, steps_tmp = 0;
    Matrix<double, Eigen::Dynamic, 3> velocities;
    double kinE = 0.0, kb = 1.0/11604.5;
    double kinT = 0.0, sumT = 0.0, sumT2 = 0.0, avgT, varT;

    Matter *mdbuff[check_steps];
    for(long i =0; i < check_steps;i++){
    mdbuff[i] = new Matter(parameters);
    }

    Dynamics PRdynamics(reactant,parameters);
    BondBoost Bbm(reactant,parameters);
    if(boost){   
        Bbm.initial();
    }

    PRdynamics.velocityScale(parameters->mdTemperature);
    dephase();

    printf("\nStarting MD run\nTemperature: %.2f Kelvin\nTotal Time: %.2f fs\nTime Step: %.2f fs\n\n",
    parameters->mdTemperature, 10.18*parameters->mdSteps*parameters->mdTimeStep, 10.18*parameters->mdTimeStep);

    long tenthSteps = parameters->mdSteps/10;
    //This prevents and edge case division by zero if mdSteps is < 10
    if (tenthSteps == 0) {
        tenthSteps = parameters->mdSteps;
    }
    while(!stoped){

        if(boost && !newstate){
            SPtime += Bbm.boost();
        }
        else{ SPtime += 10.18*parameters->mdTimeStep;}
                
        kinE = reactant->getKineticEnergy();
        kinT = (2*kinE/nFreeCoord/kb); 
        sumT += kinT;
        sumT2 += kinT*kinT;

        PRdynamics.oneStep(parameters->mdTemperature);

        md_fcalls++;
        ncheck++;
        nsteps++;

        if(parameters->mdRefine && remember && !newstate ){
            *mdbuff[ncheck-1] = *reactant;
          //  stepsbuff[ncheck-1] = nsteps;
            SPtimebuff[ncheck-1] = SPtime;
        }

        //printf("MDsteps %ld Ekin = %lf Tkin = %lf \n",nsteps,EKin,TKin); 

#ifndef NDEBUG
        if (ncheck == check_steps && !newstate){
           reactant->matter2xyz("movie", true);
        }
#endif
        if (ncheck == check_steps && !newstate){
            ncheck = 0; // reinitialize the ncheck
            status = checkState(reactant);
            if(status == true){
                remember = false;
            }
        }

        if (status && !newstate){
            nexam++;
            if (nexam >= relax_steps){
                nexam = 0;
                ncheck = 0; 	
                newstate = checkState(reactant);
                //stoped = newstate;
                status = false;
                if(newstate == false){
                    remember = true;
                }else{
                    *transition = *reactant;
                    steps_tmp = nsteps;
                    //printf("steps_tmp = %ld\n",steps_tmp);
                    //nsteps_refined = nsteps + 1;
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

    //Here we use Binary Search to refine the result; 	
    //for(long i =0; i < check_steps;i++){
    //  printf("%ld refine steps %ld\n",i,stepsbuff[i]);
    //}
    //nsteps = nsteps + 1;

    if(parameters->mdRefine && newstate){

        nsteps_refined = Refine(mdbuff,check_steps);
        printf("nsteps_refined=%ld\n",nsteps_refined); 
        long final_refined = steps_tmp-check_steps-relax_steps+nsteps_refined+1;
        //long totsteps = nsteps-check_steps-relax_steps+nsteps_refined+1; 
        printf("final_step = %ld\n",final_refined);
        *reactant = *mdbuff[nsteps_refined-1];
        *transition = *reactant;
        SPtime = SPtimebuff[nsteps_refined-1];

        /*
        for(long i = 0; i<relax_steps;i++){
            PRdynamics.oneStep(parameters->mdTemperature);
            totsteps ++;
            RLtime += 10*parameters->mdTimeStep;
            md_fcalls ++;
        }
        */
    }
    return;

    delete transition;
    for(long i =0; i < check_steps;i++){
        delete[] mdbuff[i];
    }
}

bool ParallelReplicaJob::checkState(Matter *matter)
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


bool ParallelReplicaJob::checkState_nq(Matter *matter) // checkstate without quench
{
    double distance = 0.0, D_tmp;
    double MoveCut = 2.0;
    Matter tmp(parameters);
    tmp = *matter;
    long nAtoms = tmp.numberOfAtoms();

    for(long int i=0;i<nAtoms;i++){
        if(!tmp.getFixed(i)){
            D_tmp = tmp.distance(*min1,i);
            //printf("number of Atom = %ld, Displacement = %lf\n",i,D_tmp);
            if(D_tmp >= MoveCut){
                distance += D_tmp;
            }
        }
    }

//     printf("Total Moved Distance = %lf\n",distance);
#ifndef NDEBUG
    printf("Total Moved Distance = %lf\n",distance);
#endif
    if (distance <= parameters->mdMaxMovedDist){
        return false;
    }
    return true;
}


void ParallelReplicaJob::saveData(int status,int bundleNumber)
{
    FILE *fileResults, *fileReactant, *fileProduct, *fileSaddle;

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
    //fprintf(fileResults, "%e total_physical_time\n", (SPtime+RLtime)*1e-15);
    fprintf(fileResults, "%e transition_time\n", SPtime*1e-15);
    //fprintf(fileResults, "%e relax_time\n", RLtime*1e-15);
    fprintf(fileResults, "%ld random_seed\n", parameters->randomSeed);
    fprintf(fileResults, "%lf potential_energy_reactant\n", min1->getPotentialEnergy());
    fprintf(fileResults, "%lf potential_energy_product\n", min2->getPotentialEnergy());
    fprintf(fileResults, "%ld potential_tag\n", parameters->potentialTag);
    fprintf(fileResults, "%ld total_force_calls\n", total_fcalls);
    fprintf(fileResults, "%ld force_calls_dephase\n", dh_fcalls);
    fprintf(fileResults, "%ld force_calls_dynamics\n", md_fcalls);
    fprintf(fileResults, "%ld force_calls_minimization\n", min_fcalls);
    fprintf(fileResults, "%ld force_calls_refine\n", rf_fcalls);

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
        snprintf(filename, STRING_SIZE, "product_%i.con", bundleNumber);
    }else{
        strncpy(filename, "product.con", STRING_SIZE);
    }

    fileProduct = fopen(filename, "wb");
    min2->matter2con(fileProduct);
    fclose(fileProduct);
  
    if (bundleNumber != -1) {
        snprintf(filename, STRING_SIZE, "saddle_%i.con", bundleNumber);
    }else{
        strncpy(filename, "saddle.con", STRING_SIZE);
    }

    fileSaddle = fopen(filename, "wb");
    transition->matter2con(fileSaddle);
    fclose(fileSaddle);

    return;
}


long ParallelReplicaJob::Refine(Matter *buff[],long length)
{
    long a1, b1, test, refined , initial, final, diff, RefineAccuracy;
    long tmp_fcalls;
    bool ytest;

    RefineAccuracy = parameters->mdRefineAccuracy; 
    printf("Starting search for transition step with accuracy of %ld steps\n", RefineAccuracy);
    ytest = false;

    initial = 0;
    final = length - 1;
    a1 = initial;
    b1 = final;
    diff = final - initial;
    test = int((b1-a1)/2);
    //printf("diff = %ld , ReAcc = %ld\n", diff,RefineAccuracy);
  
    tmp_fcalls = min_fcalls ;
    min_fcalls = 0;
    while(diff > RefineAccuracy)
    {
        test = a1+int((b1-a1)/2);
        ytest = checkState(buff[test]);

        if ( ytest == 0 ){
            a1 = test;
            b1 = b1;
        }
        else if ( ytest == 1 ){
            a1 = a1;
            b1 = test;
        }
        else { 
            printf("Refine Step Failed ! \n");
            exit(1);
        }
        diff = abs( b1 - a1 );
    //   printf("Insert Point %ld; Test ytest = %d ; New Boundary [ %ld, %ld ] \n",test,ytest,a1,b1);
    }

    rf_fcalls = min_fcalls;
    min_fcalls = tmp_fcalls;

    refined = int((a1+b1)/2);
    //printf("Refined mdsteps = %ld\n",nsteps_refined);
    return refined;
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
//        printf("nbuff = %ld\n",nbuff);
        Matter *DHbuff[nbuff];
        for(i=0; i<nbuff; i++){
           DHbuff[i] = new Matter(parameters);
        }

        for(i=0; i<nbuff; i++){
           dephaseDynamics.oneStep(parameters->mdTemperature);
           dh_fcalls++;
           *DHbuff[i] = *reactant;
        }
        
        state = checkState(reactant);
//       printf("state = %d\n",state);  
        if(state == true){
            dh_refined = Refine(DHbuff,nbuff);
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
}

