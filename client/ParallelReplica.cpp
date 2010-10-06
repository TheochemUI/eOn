#include "Matter.h"
#include "Constants.h"
#include "Dynamics.h"
#include "ParallelReplica.h"
#include "ConjugateGradients.h"

ParallelReplica::ParallelReplica(Parameters *params)
{
    parameters = params;
    nsteps = 0;
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

    stoped = false;  // initial stoped to false
    dynamics();
    saveData(newstate,bundleNumber);
    
    if(newstate){
        printf("New state has been found with %ld steps !\n", nsteps);
    }else{
       printf("New state has not been found in this %ld Dynamics steps !\n",parameters->mdSteps);
    }

    delete min1;
    delete min2;
    delete reactant;
}

void ParallelReplica::dynamics()
{
	bool   status = false;
    long   nFreeCoord_, nexam = 0;
    double *freeVelocities;
    double EKin=0.0, kb = 1.0/11604.5;
    double TKin=0.0, SumT = 0.0, SumT2 = 0.0, AvgT, VarT;
       
    nFreeCoord_=3*reactant->numberOfFreeAtoms();
    freeVelocities = new double[nFreeCoord_];

    Dynamics PRdynamics(reactant,parameters);
        
    PRdynamics.velocityScale();

    while(!stoped)
    {
        PRdynamics.oneStep();
        nsteps++;
        md_fcalls++;
        
        reactant->getFreeVelocities(freeVelocities);
        EKin = reactant->kineticEnergy();
        TKin = (2*EKin/nFreeCoord_/kb); 
        SumT += TKin;
        SumT2 += TKin*TKin;
        //printf("MDsteps %ld Ekin = %lf Tkin = %lf \n",nsteps,EKin,TKin); 
        
        if (nsteps % parameters->CheckFreq == 0){
#ifndef NDEBUG           
            reactant->matter2xyz("movie", true);
#endif
			status = firstAchieve();
        }
       
		if (status){
			nexam ++;
		    if (nexam >= parameters->NewRelaxSteps){	
				newstate = IsNewState();
				nexam = 0;
			}
		}

        if (nsteps >= parameters->mdSteps ){
           stoped = true;
        }       
    }
     
    AvgT=SumT/nsteps;
    VarT=SumT2/nsteps-AvgT*AvgT;
    printf("\nTempeture : Average = %lf ; Variance = %lf ; Factor = %lf \n", AvgT,VarT,VarT/AvgT/AvgT*nFreeCoord_/2);
    

    delete [] freeVelocities;
    return;
};

bool ParallelReplica::firstAchieve()
{
	double distance; 
 	*min2 = *reactant;

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
    if (distance <= parameters->PRD_MaxMovedDist){
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
    fprintf(fileResults, "%ld random_seed\n", parameters->randomSeed);
    fprintf(fileResults, "%ld potential_tag\n", parameters->potentialTag);
    fprintf(fileResults, "%ld total_force_calls\n", total_fcalls);
    fprintf(fileResults, "%ld force_calls_minimization\n", min_fcalls);
    fprintf(fileResults, "%ld force_calls_dynamics\n", md_fcalls);
    fprintf(fileResults, "%lf moved_distance\n",
            min2->distanceTo(*min1));
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


