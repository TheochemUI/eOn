//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//
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
    check_steps = parameters->CheckFreq;
    relax_steps = parameters->NewRelaxSteps;
    SPtimebuff = new double[check_steps];
    newstate = false;
    min_fcalls = 0;
    md_fcalls = 0;
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
    cgMin1.fullRelax();
    min_fcalls += min1->getForceCalls();
        
    printf("Now running Parallel Replica Dynamics\n\n");

    dynamics();

    *min2 = *reactant;   
    ConjugateGradients cgMin2(min2, parameters);
    cgMin2.fullRelax();
    min_fcalls += min2->getForceCalls();
 
    saveData(newstate,bundleNumber);
    
    if(newstate){
        printf("New state has been found\n");
        printf("Transition Time: %.2e\n", SPtime*1e-15);
    }else{
       printf("New state has not been found in this %ld dynamics steps (%.2f fs)\n",
            parameters->mdSteps,10*parameters->mdSteps*parameters->mdTimeStep);
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
    bool   boost = parameters->BondBoost;
    long   nFreeCoord = reactant->numberOfFreeAtoms()*3;
    long   ncheck = 0, nexam = 0, steps_tmp = 0;
    Matrix<double, Eigen::Dynamic, 3> velocities;
    double EKin=0.0, kb = 1.0/11604.5;
    double TKin=0.0, SumT = 0.0, SumT2 = 0.0, AvgT, VarT;
    
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
           parameters->mdTemperature, 10*parameters->mdSteps*parameters->mdTimeStep,10*parameters->mdTimeStep);
    while(!stoped){

        if(boost && !newstate){
           SPtime += Bbm.boost();
        }
        else{ SPtime += 10*parameters->mdTimeStep;}
                
        EKin = reactant->getKineticEnergy();
        TKin = (2*EKin/nFreeCoord/kb); 
        SumT += TKin;
        SumT2 += TKin*TKin;

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
  	    ncheck = 0; // reinitial the ncheck
            status = CheckState(reactant);
            if(status == true){
               remember = false;
            }
        }
       
        if (status && !newstate){
                nexam ++;
            if (nexam >= relax_steps){
                    nexam = 0;
                    ncheck = 0; 	
                    newstate = CheckState(reactant);
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

  		
        if (nsteps >= parameters->mdSteps ){
           stoped = true;
        }

        if (nsteps % 500 == 0) {
            // Since we only have a bundle size of 1 we can play with boinc_fraction_done
            // directly. When we have done parameters->mdSteps number of steps we aren't 
            // quite done so I increase the max steps by 5%.
            boinc_fraction_done((double)nsteps/(double)(parameters->mdSteps+0.05*parameters->mdSteps));
        }

        if (nsteps % (parameters->mdSteps/10) == 0 || nsteps == parameters->mdSteps) {
            printf("progress: %3.0f%%, step %7ld/%ld\n", (double)100.0*nsteps/parameters->mdSteps,
                   nsteps, parameters->mdSteps);
        }

    }
    AvgT=SumT/nsteps;
    VarT=SumT2/nsteps-AvgT*AvgT;
    printf("\nTemperature : Average = %lf ; Variance = %lf ; Factor = %lf\n\n", AvgT,VarT,VarT/AvgT/AvgT*nFreeCoord/2);

    //Here we use Binary Search to refine the result; 	
    //for(long i =0; i < check_steps;i++){
	//	printf("%ld refine steps %ld\n",i,stepsbuff[i]);		
	//}
    //nsteps = nsteps + 1;
    
   
    if(parameters->mdRefine && newstate){     
        
        Refine(mdbuff);
         printf("nsteps_refined=%ld\n",nsteps_refined); 
        long final_refined = steps_tmp-check_steps-relax_steps+nsteps_refined+1;
        long totsteps = nsteps-check_steps-relax_steps+nsteps_refined+1; 
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
};

bool ParallelReplicaJob::CheckState(Matter *matter)
{
     double distance; 
     Matter tmp(parameters);
     tmp = *matter;
     ConjugateGradients cgMin(&tmp, parameters);
     cgMin.fullRelax();
     min_fcalls += tmp.getForceCalls();


     distance = tmp.distanceTo(*min1);
#ifndef NDEBUG
     printf("Total Moved Distance = %lf\n",distance);
#endif
     if (distance <= parameters->PRD_MaxMovedDist){
        return false;
     };
     return true;
}

bool ParallelReplicaJob::CheckState_nq(Matter *matter) // checkstate without quench
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
     if (distance <= parameters->PRD_MaxMovedDist){
        return false;
     };
     return true;
}




void ParallelReplicaJob::saveData(int status,int bundleNumber){
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
     long total_fcalls = min_fcalls + md_fcalls;

     fprintf(fileResults, "%d termination_reason\n", status);
     //fprintf(fileResults, "%e total_physical_time\n", (SPtime+RLtime)*1e-15);
     fprintf(fileResults, "%e transition_time\n", SPtime*1e-15);
     //fprintf(fileResults, "%e relax_time\n", RLtime*1e-15);
     fprintf(fileResults, "%ld random_seed\n", parameters->randomSeed);
     fprintf(fileResults, "%lf potential_energy_reactant\n", min1->getPotentialEnergy());
     fprintf(fileResults, "%lf potential_energy_product\n", min2->getPotentialEnergy());
     fprintf(fileResults, "%ld potential_tag\n", parameters->potentialTag);
     fprintf(fileResults, "%ld total_force_calls\n", total_fcalls);
     fprintf(fileResults, "%ld force_calls_minimization\n", min_fcalls);
     fprintf(fileResults, "%ld force_calls_dynamics\n", md_fcalls);
     fprintf(fileResults, "%lf moved_distance\n",min2->distanceTo(*min1));
     fclose(fileResults);

     if (bundleNumber != -1) {
         snprintf(filename, STRING_SIZE, "reactant_%i.con", bundleNumber);
     }else{
         strncpy(filename, "reactant.con", STRING_SIZE);
     }
     fileReactant = fopen(filename, "wb");
     min1->matter2con(fileReactant);

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


void ParallelReplicaJob::Refine(Matter *mdbuff[]){
   
     long a1, b1, test , initial, final, diff, RefineAccuracy;
     bool ytest;

     RefineAccuracy = parameters->RefineAccuracy; 
     printf("Starting search for transition step with accuracy of %ld steps\n", RefineAccuracy);
     ytest = false;
     
     initial = 0;
     final = check_steps - 1;
     a1 = initial;
     b1 = final;
     diff = final - initial;
     test = int((b1-a1)/2);
     //printf("diff = %ld , ReAcc = %ld\n", diff,RefineAccuracy);

     while(diff > RefineAccuracy){

         test = a1+int((b1-a1)/2);
         ytest = CheckState(mdbuff[test]);
         
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
     //    printf("Insert Point %ld; Test ytest = %d ; New Bondary [ %ld, %ld ] \n",test,ytest,a1,b1);
     }

    
     nsteps_refined = int((a1+b1)/2);
     //printf("Refined mdsteps = %ld\n",nsteps_refined);
     return;
}

void ParallelReplicaJob::dephase(){

     long  DHsteps = parameters->DephaseSteps,i = 0,step_buff = 0;
     long  scType = parameters->DH_CheckType;
     long  DHcorrect = parameters->DephaseConstrain;
     bool  state = false, stop = false;
     Matrix<double, Eigen::Dynamic, 3> velocity;      
     Matter *initial;
     Matter *buff;
     initial = new Matter(parameters);      
     buff = new Matter(parameters);
     *initial = *reactant;
     *buff = *reactant;

//     printf("scType = %ld\n",scType);
//    printf("DH constrain = %ld\n", DHcorrect);
    
     Dynamics DHdynamics(reactant,parameters);
     printf("Dephasing for %ld steps\n",DHsteps);

     while(!stop){             

           DHdynamics.oneStep(parameters->mdTemperature);
           md_fcalls++;
          
           if(i % 1000 == 0){

              if(scType == 1){
     	         state = CheckState(reactant);
              }else if(scType == 2){
                 state = CheckState_nq(reactant);
              }else { 
  	         printf("Unknown CheckState method in Dephase Step, Use default value\n");
                 state = CheckState(reactant);
	      }

//              reactant->matter2xyz("movie", true);

       	      if(state == true){
                 if(DHcorrect == 1){
                    printf("Dephasing Warning: Get to the new state at step %ld, Dephase again\n",i);            
                    i = 0;
                    *reactant = *initial;
                    state = false;
                 }else if(DHcorrect ==2){
                    printf("Dephasing Warning: Get to the new state at step %ld, Now inverse the Momentum and restart from step %ld\n",i,step_buff);
                    i = step_buff;
                    *reactant = *buff;
                    velocity = reactant->getVelocities();
                    velocity = velocity*(-1);
                    reactant->setVelocities(velocity);
                    state = false;
                 }else{
		       printf("Unknown Constrain method in Dephase Step, Use default value\n");
                       DHcorrect = 1;
                 }
     	      }else{
		  *buff = *reactant;
                  step_buff = i;
              }

         //     printf(" Steps = %ld ; State = %d \n", i, state);
           }

           if(i == DHsteps){  
              stop = true;
              state=CheckState(reactant);
              if(state == false){
                 printf("Dephasing successful\n");
              }else{
                 printf("Warning:Dephasing Failed, Now in a new State!\n");
              }
           }
           i++;
     }
}
