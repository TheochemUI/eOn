#include <cstdlib>

#include "Matter.h"
#include "Dynamics.h"
#include "BondBoost.h"
#include "ParallelReplicaJob.h"
#include "ConjugateGradients.h"

ParallelReplicaJob::ParallelReplicaJob(Parameters *params)
{
    parameters = params;
    nsteps = 0;
    nsteps_refined = 0;
    SPtime = 0.0;
    RLtime = 0.0;
    check_steps = parameters->CheckFreq;
    relax_steps = parameters->NewRelaxSteps;
    stepsbuff = new long[check_steps];
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

    reactant->con2matter(reactant_passed);
    *min1 = *reactant;
    *min2 = *reactant;

    ConjugateGradients cgMin1(min1, parameters);
    cgMin1.fullRelax();
    min_fcalls += min1->getForceCalls();
        
    printf("Now running Parralel Replica Dynamics\n");

    dynamics();
    saveData(newstate,bundleNumber);
    
    printf("Total Simulated Physical Time = %lf\n",SPtime+RLtime);
    printf("Physical Thansition Time = %lf\n",SPtime);
    if(newstate){
        printf("New state has been found with %ld steps (%lf fs)!\n", nsteps_refined,SPtime);
    }else{
       printf("New state has not been found in this %ld Dynamics steps (%lf fs) !\n",parameters->mdSteps,10*parameters->mdSteps*parameters->mdTimeStep);
    }

    delete min1;
    delete min2;
    delete reactant;
}

void ParallelReplicaJob::dynamics()
{
    bool   status = false, remember = true, stoped = false;
    bool   boost = parameters->BondBoost;
    long   nFreeCoord = reactant->numberOfFreeAtoms()*3;
    long   ncheck = 0, nexam = 0;
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
        
    PRdynamics.velocityScale();

    while(!stoped){
  		
        if(boost){
           SPtime += Bbm.boost();
        }
        else{ SPtime += 10*parameters->mdTimeStep;}
        
        PRdynamics.oneStep();
        
        velocities = reactant->getVelocities();
        EKin = reactant->getKineticEnergy();
        TKin = (2*EKin/nFreeCoord/kb); 
        SumT += TKin;
        SumT2 += TKin*TKin;

        md_fcalls++;
        ncheck++;
        nsteps++;

	if(parameters->mdRefine && remember){
            *mdbuff[ncheck-1] = *reactant;
            stepsbuff[ncheck-1] = nsteps;
            SPtimebuff[ncheck-1] = SPtime;
	}

        //printf("MDsteps %ld Ekin = %lf Tkin = %lf \n",nsteps,EKin,TKin); 
        

#ifndef NDEBUG           
        if (ncheck == check_steps){            
           reactant->matter2xyz("movie", true);
        }
#endif
        if (ncheck == check_steps){
  	    ncheck = 0; // reinitial the ncheck
            status = CheckState(reactant);
            if(status == true){
               remember = false;
            }
        }
       
	if (status){
            nexam ++;
	    if (nexam >= relax_steps){
                nexam = 0;
                ncheck = 0; 	
		        newstate = CheckState(reactant);
                stoped = newstate;
                status = false;
                remember = true;
            }
	}

	//printf("%ld  ncheck = %ld nexam = %ld\n ",nsteps,ncheck,nexam);

        if (nsteps >= parameters->mdSteps ){
           stoped = true;
        }       
    }
     
    AvgT=SumT/nsteps;
    VarT=SumT2/nsteps-AvgT*AvgT;
    printf("Tempeture : Average = %lf ; Variance = %lf ; Factor = %lf \n", AvgT,VarT,VarT/AvgT/AvgT*nFreeCoord/2);

    //Here we use Binary Search to refine the result; 	
    //for(long i =0; i < check_steps;i++){
	//	printf("%ld refine steps %ld\n",i,stepsbuff[i]);		
	//}
    nsteps = nsteps + 1;
    nsteps_refined = nsteps;
    if(parameters->mdRefine && newstate){     
        
        Refine(mdbuff);

        long final_refined = nsteps_refined-nsteps+check_steps+relax_steps;
        long totsteps = nsteps_refined; 

        *reactant = *mdbuff[final_refined-1];
        SPtime = SPtimebuff[final_refined-1];

        for(long i = 0; i<relax_steps;i++){
            PRdynamics.oneStep();
            totsteps ++;
            RLtime += 10*parameters->mdTimeStep;
            md_fcalls ++;
        }
    }
  
    return;
};

bool ParallelReplicaJob::CheckState(Matter *matter)
{
     double distance; 
     *min2 = *matter;

     ConjugateGradients cgMin2(min2, parameters);
     cgMin2.fullRelax();
     min_fcalls += min2->getForceCalls();

     distance = min2->distanceTo(*min1);

#ifndef NDEBUG
     printf("Total Moved Distance = %lf\n",distance);
#endif
     if (distance <= parameters->PRD_MaxMovedDist){
	return false;
     }else { 
 	return true;
     }
}


void ParallelReplicaJob::saveData(int status,int bundleNumber){
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
     long total_fcalls = min_fcalls + md_fcalls;

     fprintf(fileResults, "%d termination_reason\n", status);
     fprintf(fileResults, "%lf total Physical time\n", SPtime+RLtime);
     fprintf(fileResults, "%lf transition_time_fs\n", SPtime);
     fprintf(fileResults, "%lf relax_time_fs\n", RLtime);
     fprintf(fileResults, "%ld random_seed\n", parameters->randomSeed);
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
     reactant->matter2con(fileProduct);
     fclose(fileProduct);

     return;
}


void ParallelReplicaJob::Refine(Matter *mdbuff[]){
   
     long a1, b1, test , initial, final, diff, RefineAccuracy;
     bool ytest;

     printf("Now started to refine the Final Point!\n");
     RefineAccuracy = parameters->RefineAccuracy; 
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
         printf("Insert Point %ld; Test ytest = %d ; New Bond [ %ld, %ld ] \n",test,ytest,a1,b1);
     }

     nsteps_refined = nsteps-check_steps-relax_steps+int((a1+b1)/2);
     printf("Refined mdsteps = %ld\n",nsteps_refined);
     return;
}

