#include <cstdlib>
#include <algorithm>

#include "Matter.h"
#include "Dynamics.h"
#include "ParallelReplicaJob.h"
#include "Log.h"
#include "HelperFunctions.h"
#include "BondBoost.h"

static const char LOG_PREFIX[] = "[ParallelReplica]";

ParallelReplicaJob::ParallelReplicaJob(Parameters *params)
{
    parameters = params;
}

ParallelReplicaJob::~ParallelReplicaJob()
{
}

std::vector<std::string> ParallelReplicaJob::run(void)
{
    //load pos.con
    reactant = new Matter(parameters);
    reactant->con2matter(helper_functions::getRelevantFile(parameters->conFilename));

    //minimize the initial reactant
    log("%s Minimizing initial position\n", LOG_PREFIX);
    reactant->relax();
    reactant->matter2con("reactant.con");

    //trajectory is the matter object for the current MD configuration
    Matter *trajectory = new Matter(parameters);
    *trajectory = *reactant;
    Dynamics dynamics(trajectory, parameters);
    BondBoost bondBoost(trajectory, parameters);

    if(parameters->biasPotential == Hyperdynamics::BOND_BOOST){
        bondBoost.initialize();
        trajectory->setBiasPotential(&bondBoost);
    }

    dephase(trajectory);

    //convert from simulation times to number of steps
    int stateCheckInterval = int(floor(parameters->parrepStateCheckInterval/parameters->mdTimeStep+0.5));
    int recordInterval = int(floor(parameters->parrepRecordInterval/parameters->mdTimeStep+0.5));

    std::vector<Matter*> MDSnapshots;
    std::vector<double> MDTimes;
    double transitionTime = 0;
    Matter transitionStructure(parameters);
    int refineForceCalls = 0;

    //Main MD loop
    double simulationTime = 0.0;
    if(parameters->biasPotential == Hyperdynamics::NONE ) {
        log("%s %8s %12s %10s %12s %12s %10s\n", LOG_PREFIX,
            "Step", "Time (s)", "KE", "PE", "TE", "KinT");
    }else{
        log("%s %8s %12s %10s %10s %12s %12s %10s\n", LOG_PREFIX,
            "Step", "Time (s)", "Boost", "KE", "PE", "TE", "KinT");
    }
    for (int step=1; step<=parameters->mdSteps; step++) {
        dynamics.oneStep();
        double boost = 1.0;
        if( parameters->biasPotential == Hyperdynamics::BOND_BOOST ) {
            double boostPotential = bondBoost.boost();   
            double kB = parameters->kB;
            boost = exp(boostPotential/kB/parameters->temperature);   
             
            simulationTime += parameters->mdTimeStep*boost;
        } else {
            simulationTime += parameters->mdTimeStep;
        }

        double kinE = trajectory->getKineticEnergy();
        double potE = trajectory->getPotentialEnergy();
        double kinT = (2.0*kinE/(trajectory->numberOfFreeAtoms()*3)/parameters->kB);

        if (step % parameters->writeMoviesInterval == 0) {
            if(parameters->biasPotential == Hyperdynamics::NONE ) {
                log("%s %8ld %12.4e %10.4f %12.4f %12.4f %10.2f\n", LOG_PREFIX, 
                        step, simulationTime*parameters->timeUnit*1e-15,
                        kinE, potE, kinE+potE, kinT);
            }else{
                double boostPotential = bondBoost.boost();   
                log("%s %8ld %12.4e %10.3e %10.4f %12.4f %12.4f %10.2f\n", LOG_PREFIX,
                        step, simulationTime*parameters->timeUnit*1e-15,
                        boost, kinE, potE+boostPotential, kinE+potE+boostPotential, kinT);
            }
        }

        //Snapshots of the trajectory used for the refinement
        if (step % recordInterval == 0 && parameters->parrepRefineTransition) {
            Matter *tmp = new Matter(parameters);    
            *tmp = *trajectory;
            MDSnapshots.push_back(tmp);
            MDTimes.push_back(simulationTime);
        }

        //check for a transition every stateCheckInterval or at the end of the simulation
        if (step % stateCheckInterval == 0 || step == parameters->mdSteps) {
            log("%s Checking for transition\n", LOG_PREFIX);

            Matter min(parameters);
            min = *trajectory;
            min.relax();

            //only check for a transition if one has yet to occur 
            if (!min.compare(reactant) && transitionTime == 0) {
                log("%s Transition occurred\n", LOG_PREFIX);

                //perform the binary search for the transition structure
                if (parameters->parrepRefineTransition) {
                    log("%s Refining transition time\n", LOG_PREFIX);
                    int tmpFcalls = Potential::fcalls;
                    int snapshotIndex = refineTransition(MDSnapshots);

                    refineForceCalls += Potential::fcalls - tmpFcalls;

                    transitionTime = MDTimes[snapshotIndex];
                    transitionStructure = *MDSnapshots[snapshotIndex];

                //if not using refinement use the current configuration as the
                //transition structure
                }else{
                    transitionStructure = *trajectory;
                    transitionTime = simulationTime;
                }
                log("%s Transition time: %.3e s\n", LOG_PREFIX, transitionTime*parameters->timeUnit*1e-15);

            //at the end of the simulation perform the refinement if it hasn't happened yet
            //this ensures that if a transition isn't seen that the same number of force
            //calls will be performed on average
            }else if (step+1 == parameters->mdSteps && transitionTime == 0) {

                //fake refinement
                if (parameters->parrepRefineTransition) {
                    log("%s Simulation ended without seeing a transition\n", LOG_PREFIX);
                    log("%s Refining anyways to prevent bias...\n", LOG_PREFIX);
                    int tmpFcalls = Potential::fcalls;
                    refineTransition(MDSnapshots, true);
                    refineForceCalls += Potential::fcalls - tmpFcalls;
                }
                transitionStructure = *trajectory;
            }

            for (unsigned int i=0;i<MDSnapshots.size();i++) {
                delete MDSnapshots[i];
            }
            MDSnapshots.clear();
            MDTimes.clear();
        }
    }

    //start the decorrelation dynamics from the transition structure
    int decorrelationSteps = int(floor(parameters->parrepCorrTime/parameters->mdTimeStep+0.5));
    log("%s Decorrelating: %i steps\n", LOG_PREFIX, decorrelationSteps);
    for (int step=1; step<=decorrelationSteps; step++) {
        dynamics.oneStep(step);
    }
    log("%s Decorrelation complete\n", LOG_PREFIX);

    //minimize the final structure
    Matter product(parameters);
    product = *trajectory;
    product.relax();
    product.matter2con("product.con");

    //report the results
    FILE *fileResults;
    std::string resultsFilename("results.dat");
    returnFiles.push_back(resultsFilename);
    fileResults = fopen(resultsFilename.c_str(), "wb");
    fprintf(fileResults, "%s potential_type\n", parameters->potential.c_str());
    fprintf(fileResults, "%ld random_seed\n", parameters->randomSeed);
    fprintf(fileResults, "%f potential_energy_reactant\n", reactant->getPotentialEnergy());
    fprintf(fileResults, "%i force_calls_refine\n", refineForceCalls);
    fprintf(fileResults, "%d total_force_calls\n", Potential::fcalls);

    if (transitionTime == 0) {
        fprintf(fileResults, "0 transition_found\n");
        fprintf(fileResults, "%e simulation_time_s\n", simulationTime*parameters->timeUnit*1.0e-15);
    }else{
        fprintf(fileResults, "1 transition_found\n");
        fprintf(fileResults, "%e transition_time_s\n", transitionTime*parameters->timeUnit*1.0e-15);
        fprintf(fileResults, "%e correlation_time_s\n", parameters->parrepCorrTime*parameters->timeUnit*1.0e-15);
        fprintf(fileResults, "%lf potential_energy_product\n", product.getPotentialEnergy());
    }
    fprintf(fileResults, "%lf speedup\n", simulationTime/(parameters->mdSteps*parameters->mdTimeStep));

    fclose(fileResults);

    MDSnapshots.clear();
    MDTimes.clear();
    delete trajectory;
    delete reactant;

    return returnFiles;
}

void ParallelReplicaJob::dephase(Matter *trajectory)
{
    Dynamics dynamics(trajectory, parameters);

    int dephaseSteps = int(floor(parameters->parrepDephaseTime/parameters->mdTimeStep+0.5));
    log("%s Dephasing: %i steps\n", LOG_PREFIX, dephaseSteps);

    Matter initial(parameters);
    initial = *trajectory;

    while (true) {
        //always start from the initial configuration
        *trajectory = initial;
        dynamics.setThermalVelocity();

        // Dephase MD trajectory
        for (int step=1; step<=dephaseSteps; step++) {
            dynamics.oneStep(step);
        }

        // Check to see if a transition occured
        Matter min(parameters);
        min = *trajectory;
        min.relax();

        if (min.compare(reactant)) {
            log("%s Dephasing successful\n", LOG_PREFIX);
            break;
        }else{
            log("%s Transition occured during dephasing; Restarting\n", LOG_PREFIX);
        }
    }
}

int ParallelReplicaJob::refineTransition(std::vector<Matter*> MDSnapshots, bool fake)
{
    int min, max, mid;
    bool midTest;
    min = 0;
    max = MDSnapshots.size() - 1;

    while( (max-min) > 1 ) {
        mid = min + (max-min)/2;
        Matter *snapshot = MDSnapshots[mid];
        snapshot->relax(true);

        if (fake == false) {
            midTest = snapshot->compare(reactant);
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

