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
#include "ReplicaExchangeJob.h"
#include "ConjugateGradients.h"

ReplicaExchangeJob::ReplicaExchangeJob(Parameters *params)
{
    parameters = params;
}

ReplicaExchangeJob::~ReplicaExchangeJob(){ }

void ReplicaExchangeJob::run(int bundleNumber)
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
    reactant->con2matter(reactant_passed);

    long NumofReplica,nsteps,i,j,nFreeCoord;
    double BaseTemp, EKin, TKin,kb,AvgT,VarT;

    NumofReplica = 10;
    nsteps = parameters->mdSteps;
    nFreeCoord = reactant->numberOfFreeAtoms()*3;
    kb = 1.0/11604.5;
    BaseTemp = parameters->mdTemperature;
    printf("BaseTemp = %lf \n",BaseTemp);    
    ReplicaTemp = new double[NumofReplica];
    Dynamics *RES[NumofReplica];
    double SumT[NumofReplica],SumT2[NumofReplica];


    for(i=0; i < NumofReplica; i++){
        ReplicaArray[i] = new Matter(parameters);
        *ReplicaArray[i] = *reactant;
        printf("%ld replica : numAtoms = %ld \n", i, ReplicaArray[i]->numberOfAtoms());
        ReplicaTemp[i] = BaseTemp + i*30;
        RES[i]=new Dynamics(ReplicaArray[i],parameters);
        RES[i]->velocityScale(ReplicaTemp[i]);
    }

    printf("Now running Replica Exchange Simulation\n");
    long numAtoms = 0;
    for(i=0; i<NumofReplica; i++){
        numAtoms = ReplicaArray[i]->numberOfAtoms();
        printf("%ld replica : numAtoms = %ld \n", i, numAtoms);
        printf("%ld replica : Temperature = %lf \n", i, ReplicaTemp[i]);
    }

    for(j= 0; j< nsteps; j++){
        printf("Steps %ld", j);
        for(i= 1; i< NumofReplica; i++){
            EKin = ReplicaArray[i]->getKineticEnergy();
            TKin = (2*EKin/nFreeCoord/kb);
            SumT[i] += TKin;
            SumT2[i] += TKin*TKin;
            printf(" %5.2lf ",TKin); 
            RES[i]->oneStep(ReplicaTemp[i]);
        }
       printf("\n");
    }
 
    for(i=1;i<NumofReplica;i++){   
        AvgT=SumT[i]/nsteps;
        VarT=SumT2[i]/nsteps-AvgT*AvgT;
        printf("Replica %ld -- %5.2lf: AvgT = %lf ; VarT = %lf ; Factor = %lf \n",i, ReplicaTemp[i],AvgT,VarT,VarT/AvgT/AvgT*nFreeCoord/2);
    }

    printf("%ld replica : Temperature = %lf \n", i, ReplicaTemp[2]);
    RES[2]->fullSteps(ReplicaTemp[2]);

    delete reactant;
    delete  ReplicaArray[2];
}

