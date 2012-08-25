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
#include <algorithm>

#include "Matter.h"
#include "Dynamics.h"
#include "UnbiasedParallelReplicaJob.h"
#include "Log.h"
#include "HelperFunctions.h"

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

static const char LOG_PREFIX[] = "[UnbiasedParallelReplica]";

UnbiasedParallelReplicaJob::UnbiasedParallelReplicaJob(Parameters *parameters_passed)
{
    parameters = parameters_passed;
}

UnbiasedParallelReplicaJob::~UnbiasedParallelReplicaJob()
{
    delete reactant;
}

std::vector<std::string> UnbiasedParallelReplicaJob::run(void)
{
    reactant = new Matter(parameters);
    reactant->con2matter(helper_functions::getRelevantFile(parameters->conFilename));

    //Always minimize the initial reactant
    log("%s Minimizing initial reactant\n", LOG_PREFIX);
    reactant->relax();
    reactant->matter2con("reactant.con");

    //trajectory is the matter object for the current MD configuration
    Matter *trajectory = new Matter(parameters);
    *trajectory = *reactant;

    dephase(trajectory);

    //convert from simulation times to number of steps
    int stateCheckInterval = int(parameters->parrepStateCheckInterval/parameters->mdTimeStepInput);
    int recordInterval = int(parameters->parrepRecordInterval/parameters->mdTimeStepInput);

    Dynamics dynamics(trajectory, parameters);
    std::vector<Matter*> MDSnapshots;
    std::vector<double> MDTimes;
    double transitionTime = 0;
    Matter transitionStructure(parameters);

    //Main MD loop
    for (int step=1;step<=parameters->mdSteps;step++) {
        dynamics.oneStep();
        double simulationTime = step*parameters->mdTimeStepInput;

        //Snapshots of the trajectory used for the refinement
        if (step % recordInterval == 0 && parameters->parrepRefineTransition) {
            Matter *tmp = new Matter(parameters);    
            *tmp = *trajectory;
            MDSnapshots.push_back(tmp);
            MDTimes.push_back(simulationTime);
        }

        //check for transition every stateCheckInterval or at the end of the simulation
        if (step % stateCheckInterval == 0 || step == parameters->mdSteps) {
            log("%s checking for transition\n", LOG_PREFIX);

            Matter min(parameters);
            min = *trajectory;
            min.relax(true);

            //only check for a transition if one has yet to occur 
            if (min != *reactant && transitionTime == 0) {
                log("%s transition occurred\n", LOG_PREFIX);
                log("%s refining transition time\n", LOG_PREFIX);

                //perform the binary search for the transition structure
                if (parameters->parrepRefineTransition) {
                    int snapshotIndex = refineTransition(MDSnapshots);
                    transitionTime = MDTimes[snapshotIndex];
                    transitionStructure = *MDSnapshots[snapshotIndex];

                //use the current configuration as the transition structure
                }else{
                    transitionStructure = *trajectory;
                    transitionTime = simulationTime;
                }
                log("%s transition occurred at %.3f\n", LOG_PREFIX, transitionTime);

            //at the end of the simulation perform the refinement if it hasn't happened yet
            //this ensures that if a transition isn't seen that the same number of force
            //calls will be performed on average
            }else if (step == parameters->mdSteps && transitionTime == 0) {

                //fake refinement
                if (parameters->parrepRefineTransition) {
                    log("%s simulation ended without seeing a transition\n", LOG_PREFIX);
                    log("%s refining anyways to prevent bias...\n", LOG_PREFIX);
                    refineTransition(MDSnapshots, true);
                }
                transitionStructure = *trajectory;

            }

            MDSnapshots.clear();
            MDTimes.clear();
        }

    }

    //start the decorrelation dynamics from the transition structure
    *trajectory = transitionStructure;
    int decorrelationSteps = int(parameters->parrepCorrTime/parameters->mdTimeStepInput);
    log("%s decorrelating for %i steps\n", LOG_PREFIX, decorrelationSteps);
    for (int step=1;step<=decorrelationSteps;step++) {
        dynamics.oneStep();
    }

    //minimize the final structure
    trajectory->relax(true);
    trajectory->matter2con("product.con");

    FILE *fileResults;
    std::string resultsFilename("results.dat");
    returnFiles.push_back(resultsFilename);
    fileResults = fopen(resultsFilename.c_str(), "wb");
    fprintf(fileResults, "%s potential_type\n", parameters->potential.c_str());
    fprintf(fileResults, "%ld random_seed\n", parameters->randomSeed);
    fprintf(fileResults, "%lf potential_energy_reactant\n", reactant->getPotentialEnergy());
    fprintf(fileResults, "%d total_force_calls\n", Potential::fcalls);


    if (transitionTime == 0) {
        fprintf(fileResults, "0 transition_found\n");
        fprintf(fileResults, "%e simulation_time_s\n", parameters->mdTime*1.0e-15);
    }else{
        fprintf(fileResults, "%e transition_time_s\n", transitionTime*1.0e-15);
        fprintf(fileResults, "%e corr_time_s\n", parameters->parrepCorrTime*1.0e-15);
        fprintf(fileResults, "%lf potential_energy_product\n", trajectory->getPotentialEnergy());
    }


    MDSnapshots.clear();
    MDTimes.clear();
    delete trajectory;

    return returnFiles;
}

void UnbiasedParallelReplicaJob::dephase(Matter *trajectory)
{
    Dynamics dynamics(trajectory, parameters);
    dynamics.setThermalVelocity();

    int dephaseSteps = int(parameters->parrepDephaseTime/parameters->mdTimeStepInput);
    log("%s dephasing for %i steps\n", LOG_PREFIX, dephaseSteps);

    while (true) {
        // Dephase MD trajectory
        for (int step=1;step<=dephaseSteps;step++) {
            dynamics.oneStep();
        }

        // Check to see if a transition occured
        Matter min(parameters);
        min = *trajectory;
        min.relax(true);
        if (min == *reactant) {
            // Dephasing successful
            log("%s dephasing successful\n", LOG_PREFIX);
            break;
        }else{
            // A transition occurred. Retry.
            log("%s transition occured during dephasing, restarting\n", LOG_PREFIX);
            continue;
        }
    }
}

int UnbiasedParallelReplicaJob::refineTransition(std::vector<Matter*> MDSnapshots, bool fake)
{
    long min, max, mid;
    bool midTest;
    min = 0;
    max = MDSnapshots.size() - 1;

    while( (max-min) > 1 ) {

        mid = min + (max-min)/2;
        Matter *snapshot = MDSnapshots[mid];
        snapshot->relax(true);

        if (fake == false) {
            midTest = *snapshot == *reactant;
        }else{
            //if we are faking the refinement just generate a random answer
            //for the comparison test
            midTest = bool(helper_functions::randomInt(0,1));
        }

        if (midTest){
            min = mid;
        } else {
            max = mid;
        }
    }

    return (min+max)/2 + 1;
}

