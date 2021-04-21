
#include <cstdlib>

#include "Matter.h"
#include "Log.h"
#include "Dynamics.h"
#include "ReplicaExchangeJob.h"
#include "HelperFunctions.h"

ReplicaExchangeJob::ReplicaExchangeJob(Parameters *params)
{
    parameters = params;
}

ReplicaExchangeJob::~ReplicaExchangeJob()
{
}

std::vector<std::string> ReplicaExchangeJob::run(void)
{
    long i, step, samplingSteps = long(parameters->repexcSamplingTime/parameters->mdTimeStep+0.5);
    long exchangePeriodSteps = long(parameters->repexcExchangePeriod/parameters->mdTimeStep+0.5);
    double energyLow, energyHigh;
    double kbTLow, kbTHigh;
    double kB = parameters->kB;
    double pAcc;
    Matter *tmpMatter;

    string posFilename = helper_functions::getRelevantFile(parameters->conFilename);
    pos = new Matter(parameters);
    pos->con2matter(posFilename);

    log("\nRunning Replica Exchange\n\n");

    long refForceCalls = Potential::fcalls;

    // allocate a Matter and Dynamics object for each replica
    Matter *replica[parameters->repexcReplicas];
    Dynamics *replicaDynamics[parameters->repexcReplicas];
    for(i=0; i<parameters->repexcReplicas; i++) {
        replica[i] = new Matter(parameters);
        *replica[i] = *pos;
        replicaDynamics[i] = new Dynamics(replica[i], parameters);
    }

    // assign temperatures
    double replicaTemperature[parameters->repexcReplicas];

    log("Temperature distribution:\n");
    if(parameters->repexcTemperatureDistribution == "linear") {
        for(i=0; i<parameters->repexcReplicas; i++)
        {
            replicaTemperature[i] = parameters->repexcTemperatureLow + double(i)/double(parameters->repexcReplicas-1.0)
                *(parameters->repexcTemperatureHigh - parameters->repexcTemperatureLow);
            replicaDynamics[i]->setTemperature(replicaTemperature[i]);
        }
    } else if(parameters->repexcTemperatureDistribution == "exponential") {
        double kTemp = log(parameters->repexcTemperatureHigh / parameters->repexcTemperatureLow)/(parameters->repexcReplicas-1.0);
        // cout <<"kTemp: "<<kTemp<<endl;
        for(i=0; i<parameters->repexcReplicas; i++)
        {
            replicaTemperature[i] = parameters->repexcTemperatureLow*exp(kTemp*i);
            replicaDynamics[i]->setTemperature(replicaTemperature[i]);
            log("replica: %ld temperature %.0f \n",i+1, replicaTemperature[i]);
        }
    }

    log("\nReplica Exchange sampling for %.0f fs; %ld steps; %ld replicas.\n",
         parameters->repexcSamplingTime*10.18, samplingSteps, parameters->repexcReplicas);

    for(step=1; step<=samplingSteps; step++)
    {
        for(i=0; i<parameters->repexcReplicas; i++)
        {
            replicaDynamics[i]->oneStep();
            cout <<"step: "<<step<<" i "<<i<<" energy: "<<replica[i]->getPotentialEnergy()<<endl;
        }
        if( (step % exchangePeriodSteps) == 0 )
        {
            for(long trial=0; trial<parameters->repexcExchangeTrials; trial++)
            {
                i = helper_functions::randomInt(0, parameters->repexcReplicas-2);
                energyLow = replica[i]->getPotentialEnergy();
                energyHigh = replica[i+1]->getPotentialEnergy();
                kbTLow = kB*replicaTemperature[i];
                kbTHigh = kB*replicaTemperature[i+1];
                pAcc = min(1.0, exp((energyHigh-energyLow)*(1.0/kbTHigh-1.0/kbTLow)));
                double tmp = helper_functions::randomDouble();
                cout <<"step: "<<step<<" trial swap, i "<<i<<" elow: "<<energyLow<<" ehigh: "<<energyHigh<<" pAcc: "<<pAcc<<" rand: "<<tmp;
                //if(helper_functions::randomDouble()<pAcc)
                if(tmp<pAcc)
                {
                    cout <<" swap\n";
                    // swap configurations
                    tmpMatter = replica[i];
                    replica[i] = replica[i+1];
                    replica[i+1] = tmpMatter;
                    // reset velocities
                    replicaDynamics[i]->setThermalVelocity();
                    replicaDynamics[i+1]->setThermalVelocity();
                }
                else { cout <<" no_swap\n"; }
            }            
        }
    }

    forceCalls = Potential::fcalls - refForceCalls;
    saveData();

    // delete Matter and Dynamics objects

    delete pos;
    return returnFiles;
}

void ReplicaExchangeJob::saveData(void)
{

    FILE *fileResults, *filePos;

    std::string resultsFilename("results.dat");
    returnFiles.push_back(resultsFilename);
    fileResults = fopen(resultsFilename.c_str(), "wb");

    fprintf(fileResults, "%ld random_seed\n", parameters->randomSeed);
    fprintf(fileResults, "%s potential_type\n", parameters->potential.c_str());
    fprintf(fileResults, "%ld force_calls_sampling\n", forceCalls);
    fclose(fileResults);

    std::string posFilename("pos_out.con");
    returnFiles.push_back(posFilename);
    filePos = fopen(posFilename.c_str(), "wb");
    fclose(filePos);

    return;
}
