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
#include "Log.h"
#include "Dynamics.h"
#include "ReplicaExchangeJob.h"
#include "HelperFunctions.h"

ReplicaExchangeJob::ReplicaExchangeJob(Parameters *parameters_passed)
{
    parameters = parameters_passed;
}

ReplicaExchangeJob::~ReplicaExchangeJob()
{
}

std::vector<std::string> ReplicaExchangeJob::run(void)
{
    long i, step, samplingSteps = long(parameters->repexcSamplingTime/parameters->mdTimeStep);
    long exchangePeriodSteps = long(parameters->repexcExchangePeriod/parameters->mdTimeStep);
    double energyLow, energyHigh;
    double kbTLow, kbTHigh, kb = 1.0/11604.5;
    double pAcc;
    Matter *tmpMatter;

    string reactantPassedFilename = helper_functions::getRelevantFile(parameters->conFilename);
    reactant = new Matter(parameters);
    reactant->con2matter(reactantPassedFilename);

    log("\nRunning Replica Exchange\n\n");

    long refForceCalls = Potential::fcalls;

    // allocate a Matter and Dynamics object for each replica
    Matter *replica[parameters->repexcReplicas];
    Dynamics *replicaDynamics[parameters->repexcReplicas];
    for(i=0; i<parameters->repexcReplicas; i++) {
        replica[i] = new Matter(parameters);
        *replica[i] = *reactant;
        replicaDynamics[i] = new Dynamics(replica[i], parameters);
    }

    // assign temperatures
    double replicaTemperature[parameters->repexcReplicas];

    cout <<"parameters->repexcTemperatureLow: "<<parameters->repexcTemperatureLow<<endl;
    cout <<"parameters->repexcTemperatureHigh: "<<parameters->repexcTemperatureHigh<<endl;
    for(i=0; i<parameters->repexcReplicas; i++)
    {
        cout <<parameters->repexcTemperatureDistribution<<endl;
        if(parameters->repexcTemperatureDistribution == "linear")
        {
            replicaTemperature[i] = parameters->repexcTemperatureLow 
              + double(i)/double(parameters->repexcReplicas - 1)
                *(parameters->repexcTemperatureHigh - parameters->repexcTemperatureLow);
        }
        else if(parameters->repexcTemperatureDistribution == "exponential")
        {
            double kTemp = log(parameters->repexcTemperatureHigh / parameters->repexcTemperatureLow)/(parameters->repexcReplicas-1);
            cout <<"kTemp: "<<kTemp<<endl;
            replicaTemperature[i] = parameters->repexcTemperatureLow*exp(kTemp*i);
        }
        cout <<"replica: "<<i+1<<" temperature: "<<replicaTemperature[i]<<endl;
        replicaDynamics[i]->setTemperature(replicaTemperature[i]);
    }

    log("Replica Exchange sampling for %f fs; %ld steps; %ld replicas.\n",
         parameters->repexcSamplingTime, samplingSteps, parameters->repexcReplicas);

    for(step=1; step<=parameters->repexcSamplingTime; step++)
    {
        for(i=0; i<parameters->repexcReplicas; i++)
        {
            replicaDynamics[i]->oneStep();
        }
        if( (step % exchangePeriodSteps) == 0 )
        {
            for(long trial=0; trial<parameters->repexcExchangeTrials; trial++)
            {
                i = helper_functions::randomInt(0, parameters->repexcReplicas-2);
                energyLow = replica[i]->getPotentialEnergy();
                energyHigh = replica[i+1]->getPotentialEnergy();
                kbTLow = kb*replicaTemperature[i];
                kbTHigh = kb*replicaTemperature[i+1];
                pAcc = min(1.0, exp((energyHigh-energyLow)/(kbTHigh-kbTLow)));
                if(helper_functions::randomDouble()<pAcc)
                {
                    // swap configurations
                    tmpMatter = replica[i];
                    replica[i] = replica[i+1];
                    replica[i+1] = tmpMatter;
                    // reset velocities
                    replicaDynamics[i]->setThermalVelocity();
                    replicaDynamics[i+1]->setThermalVelocity();
                }
            }            
        }
    }

    forceCalls = Potential::fcalls - refForceCalls;
    saveData();

    // delete Matter and Dynamics objects

    delete reactant;
    return returnFiles;
}

void ReplicaExchangeJob::saveData(void)
{

    FILE *fileResults, *fileReactant, *fileProduct;

    std::string resultsFilename("results.dat");
    returnFiles.push_back(resultsFilename);
    fileResults = fopen(resultsFilename.c_str(), "wb");

    fprintf(fileResults, "%ld random_seed\n", parameters->randomSeed);
    fprintf(fileResults, "%s potential_type\n", parameters->potential.c_str());
    fprintf(fileResults, "%ld force_calls_sampling\n", forceCalls);
    fclose(fileResults);

    std::string reactantFilename("reactant.con");
    returnFiles.push_back(reactantFilename);
    fileReactant = fopen(reactantFilename.c_str(), "wb");
    fclose(fileReactant);

    return;
}
