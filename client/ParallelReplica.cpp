#include <cstdlib>

#include "Matter.h"
#include "Constants.h"
#include "Dynamics.h"
#include "ParallelReplica.h"
#include "ConjugateGradients.h"

ParallelReplica::ParallelReplica(Parameters *params)
{
    parameters = params;
    nsteps = 0;
    ncheck = 0;
    nexam = 0;
    remember = true;
    stoped = false;
    newstate = false;
    min_fcalls = 0;
    md_fcalls = 0;
}

ParallelReplica::~ParallelReplica(){ }

void ParallelReplica::run(int bundleNumber)
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
    
    if(newstate){
        printf("New state has been found with %ld steps (%lf fs)!\n", nsteps,10*nsteps*parameters->mdTimeStep);
    }else{
       printf("New state has not been found in this %ld Dynamics steps (%lf fs) !\n",parameters->mdSteps,10*parameters->mdSteps*parameters->mdTimeStep);
    }

    delete min1;
    delete min2;
    delete reactant;
}

void ParallelReplica::dynamics()
{
    bool   status = false;
    long   nFreeCoord = reactant->numberOfFreeAtoms()*3;
    Matrix<double, Eigen::Dynamic, 3> velocities;
    double EKin=0.0, kb = 1.0/11604.5;
    double TKin=0.0, SumT = 0.0, SumT2 = 0.0, AvgT, VarT;
     
    Matter *mdbuff[parameters->CheckFreq];	
    for(long i =0; i < parameters->CheckFreq;i++){
	    mdbuff[i] = new Matter(parameters);
	}

    Dynamics PRdynamics(reactant,parameters);
        
    PRdynamics.velocityScale();

    while(!stoped){
		
        PRdynamics.oneStep();
        
        velocities = reactant->getVelocities();
        EKin = reactant->getKineticEnergy();
        TKin = (2*EKin/nFreeCoord/kb); 
        SumT += TKin;
        SumT2 += TKin*TKin;

        md_fcalls++;
        ncheck++;
        nsteps++;
	reactant->setNsteps(nsteps);


	if(parameters->mdRefine && remember){
            *mdbuff[ncheck-1] = *reactant;
	}

        //printf("MDsteps %ld Ekin = %lf Tkin = %lf \n",nsteps,EKin,TKin); 
        
        if (ncheck == parameters->CheckFreq){
#ifndef NDEBUG           
            reactant->matter2xyz("movie", true);
#endif
  	    ncheck = 0; // reinitial the ncheck
            status = firstArchieve(reactant);
        }
       
	if (status){
            nexam ++;
	    if (nexam >= parameters->NewRelaxSteps){	
		newstate = IsNewState();
            }
	}

	//printf("%ld  ncheck = %ld nexam = %ld\n ",nsteps,ncheck,nexam);

        if (nsteps >= parameters->mdSteps ){
           stoped = true;
        }       
    }
     
    AvgT=SumT/nsteps;
    VarT=SumT2/nsteps-AvgT*AvgT;
    printf("\nTempeture : Average = %lf ; Variance = %lf ; Factor = %lf \n", AvgT,VarT,VarT/AvgT/AvgT*nFreeCoord/2);

    //Here we use Golden Search to refine the result; 	
    //for(long i =0; i < parameters->CheckFreq;i++){
	//	printf("%ld refine steps %ld\n",i,mdbuff[i]->getNsteps());		
	//}
   
    if(parameters->mdRefine && newstate){     
        Refine(mdbuff);
    }

    return;
};

bool ParallelReplica::firstArchieve(Matter *matter)
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


bool ParallelReplica::IsNewState(){

     double distance;

     *min2 = *reactant;
     ConjugateGradients cgMin2(min2, parameters);
     cgMin2.fullRelax();
     min_fcalls += min2->getForceCalls();

     distance = min2->distanceTo(*min1);
#ifndef NDEBUG
     printf("Total Moved Distance = %lf\n",distance);
#endif
     nexam = 0;
     if (distance <= parameters->PRD_MaxMovedDist){
	ncheck = 0;
        return false;
     }else {
        stoped = true;
        return true;
     }
}

void ParallelReplica::saveData(int status,int bundleNumber){
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
     fprintf(fileResults, "%ld simulation_steps\n", nsteps);
     fprintf(fileResults, "%lf simulation_time_fs\n",10*nsteps*parameters->mdTimeStep);
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


void ParallelReplica::Refine(Matter *mdbuff[]){
   
     long a1, b1, a2, b2, initial, final, diff, RefineAccuracy;
     bool ya, yb;

     printf("Now started to refine the Final Point!\n");
     RefineAccuracy = parameters->RefineAccuracy; 
     ya = false;
     yb = false;
     initial = 0;
     final = parameters->CheckFreq - 1;
     a1 = a2 = initial;
     b1 = b2 = final;
     diff = final - initial;
      
     while(diff > RefineAccuracy){
	   a2 = int(a1 + 0.382 * ( b1 - a1) );
	   b2 = int(a1 + 0.618 * ( b1 - a1) );
	   ya = firstArchieve(mdbuff[a2]);
	   yb = firstArchieve(mdbuff[b2]);
     
	   if ( ya == 0 && yb == 0){
	      a1 = initial;
              b1 = a2;
	   }else if( ya == 0 && yb == 1){
              a1 = a2;
	      b1 = b2;
	   }else if( ya == 1 && yb == 1){
	      a1 = b2;
	      b1 = final;
	   }else if( ya == 1 && yb == 0){
	      printf("Warning : Recrossing happened, search ranger will defined as (b2,final)\n");
	      a1 = b2;
	      b1 = final;
	   }else {
	      exit(1);
	   }
	
	   diff = abs(b2 -a2);
     }
	
     printf("Refined mdsteps = %ld\n",nsteps-parameters->CheckFreq-parameters->NewRelaxSteps+(a2+b2)/2);
     return;
}

