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
#include "Optimizer.h"
#include "Log.h"

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
}

ParallelReplicaJob::~ParallelReplicaJob()
{
}

std::vector<std::string> ParallelReplicaJob::run(void)
{
    current = new Matter(parameters);
    reactant = new Matter(parameters);
    transition = new Matter(parameters);
    transition_relaxed = new Matter(parameters);
    product = new Matter(parameters);
    product_relaxed = new Matter(parameters);

    minimizeFCalls = mdFCalls = refineFCalls = dephaseFCalls = corrFCalls = 0;
    time = transitionTime = corrTime = 0.0;
    jobStatus = ParallelReplicaJob::STATUS_NOTRAN;
    newStateFlag = false;
    relaxStatus = true;
    string reactantFilename = helper_functions::getRelevantFile(parameters->conFilename);
    current->con2matter(reactantFilename);
    log("\nMinimizing initial reactant\n");
    long refFCalls = Potential::fcalls;
    *reactant = *current;
    relaxStatus = reactant->relax();
    if(!relaxStatus){
        jobStatus = ParallelReplicaJob::STATUS_BAD_RELAXFAILED;
    }

    minimizeFCalls += (Potential::fcalls - refFCalls);

    log("\nParallel Replica Dynamics, running\n\n");
    
    int state=dynamics();

    saveData(state);
    printEndStatus();

    if(newStateFlag){
        log("[ParallelReplica] Transition time: %.2e s\n", transitionTime*1.0e-15*parameters->timeUnit);
    }else{
       log("[ParallelReplica] No new state was found in %ld dynamics steps and a time: %.3e s\n",
           parameters->mdSteps, time*1.0e-15*parameters->timeUnit);
    }

    delete current;
    delete reactant;
    delete transition;
    delete transition_relaxed;
    delete product;
    delete product_relaxed;

    return returnFiles;
}

int ParallelReplicaJob::dynamics()
{
    bool transitionFlag = false, recordFlag = true, stopFlag = false, refineFlag = true;
    long nFreeCoord = reactant->numberOfFreeAtoms()*3;
    long mdBufferLength, refFCalls;
    long step = 0, refineStep, newStateStep = 0; // check that newStateStep is set before used
    long nCheck = 0, nCorr = 0, nRecord = 0, nboost = 0;
    long StateCheckInterval, RecordInterval, CorrSteps;
    double kinE, kinT, avgT, varT,  kb = 1.0/11604.5;
    double sumT = 0.0, sumT2 = 0.0;
    double sumboost = 0.0, boost = 0.0, boostPotential = 0.0;
    double refinedTime=0.0;
    Matter **mdBuffer;

    StateCheckInterval = int(floor(parameters->parrepStateCheckInterval/parameters->mdTimeStep+0.5));
    RecordInterval = int(floor(parameters->parrepRecordInterval/parameters->mdTimeStep+0.5));
    CorrSteps = int(floor(parameters->parrepCorrTime/parameters->mdTimeStep+0.5));
    refineFlag = parameters->parrepRefineTransition;


    mdBufferLength = long(floor(StateCheckInterval/RecordInterval+0.5));
    mdBuffer = new Matter *[mdBufferLength];
    for(long i=0; i<mdBufferLength; i++) {
        mdBuffer[i] = new Matter(parameters);
    }
    timeBuffer = new double[mdBufferLength];
  
    Dynamics parrepDynamics(current, parameters);
    BondBoost bondBoost(current, parameters);

    if(parameters->biasPotential == Hyperdynamics::BOND_BOOST){
        bondBoost.initialize();
    }

    parrepDynamics.setThermalVelocity();

    // dephase the trajectory so that it is thermal and independent of others
    refFCalls = Potential::fcalls;
    dephase();
    dephaseFCalls = Potential::fcalls - refFCalls;

    log("\nStarting MD run\nTemperature: %.2f K\nTotal simulation time: %.2f fs\nTime step: %.2f fs\nTotal steps: %ld\n\n", 
        parameters->temperature, 
        parameters->mdSteps*parameters->mdTimeStep*parameters->timeUnit,
        parameters->mdTimeStep*parameters->timeUnit,
        parameters->mdSteps);
    log("MD buffer length: %ld\n", mdBufferLength);

    long tenthSteps = parameters->mdSteps/10;
    // prevents a division by zero if mdSteps is < 10
    if (tenthSteps == 0) {
        tenthSteps = parameters->mdSteps;
    }

    // loop dynamics iterations until some condition tells us to stop
    while(!stopFlag)
    {
        if( parameters->biasPotential == Hyperdynamics::BOND_BOOST ) {
            boostPotential = bondBoost.boost();   
            boost = 1.0*exp(boostPotential/parameters->kB/parameters->temperature);   
            
            time += parameters->mdTimeStep*boost;
            if (boost > 1.0){
                sumboost += boost;
                nboost ++;
            }
            
        } else {
            time += parameters->mdTimeStep;
        }

        kinE = current->getKineticEnergy();
        kinT = (2.0*kinE/nFreeCoord/kb); 
        sumT += kinT;
        sumT2 += kinT*kinT;
        //log("steps = %10d temp = %10.5f \n",step,kinT);

        parrepDynamics.oneStep();
        mdFCalls++;

        nCheck++; // count up to parameters->parrepStateCheckInterval before checking for a transition
        step++;
        //log("step = %4d, time= %10.4f\n",step,time);
        // standard conditions; record mater object in the transition buffer
        if( parameters->parrepRefineTransition && recordFlag && !newStateFlag )
        {
            if( nCheck % RecordInterval == 0 )
            {
                *mdBuffer[nRecord] = *current;
                timeBuffer[nRecord] = time;
                nRecord++; // current location in the buffer
            }
        }

//#ifndef NDEBUG
//        if (ncheck == StateCheckInterval && !newState){
//           current->matter2xyz("movie", true);
//        }
//#endif

        // time to do a state check; if a transiton if found, stop recording
        if ( nCheck == StateCheckInterval) 
        {
            nCheck = 0; // reinitialize check state counter
            nRecord = 0; // restart the buffer
            refFCalls = Potential::fcalls;
            transitionFlag = checkState(current, reactant);
            minimizeFCalls += Potential::fcalls - refFCalls;
            // Once a transition is detected we keep on optimizing, but do nothing else;
            if(newStateFlag){
                transitionFlag = false;
            }
            //Run additional Corrstep to check for a recrossing event
            if(transitionFlag){
                transitionStep = step;
                *product = *current;
                transitionTime = timeBuffer[int(floor(mdBufferLength/2.)+0.5)];
                log("Detected transition; Running %ld MD steps to see if it will recross\n",CorrSteps);
                recordFlag = false;
            }
        }

        // a transition has been detected; run an additional relaxSteps to see if the products are stable
        if ( transitionFlag && !newStateFlag )
        {
            //if(nCorr == 0){
            //    refFCalls = Potential::fcalls;
            //}
            nCorr++; // run relaxSteps before making a state check
            if(CorrSteps > parameters->mdSteps-step){
                jobStatus = ParallelReplicaJob::STATUS_TRAN_NOTIME;
                newStateFlag = true;
                refineFlag = false;
                corrTime = nCorr*parameters->mdTimeStep;

                *product_relaxed=*product;
                relaxStatus = product_relaxed->relax(true);
                if(!relaxStatus){
                    jobStatus = ParallelReplicaJob::STATUS_BAD_RELAXFAILED;
                }
            }
            if (nCorr > CorrSteps){
                // state check; reset counters
                nCorr = 0;
                nCheck = 0;
                nRecord = 0;
                refFCalls = Potential::fcalls;
                newStateFlag = checkState(current, reactant);
                corrFCalls += Potential::fcalls - refFCalls;
                transitionFlag = false;
                if(newStateFlag == false){
                    jobStatus = ParallelReplicaJob::STATUS_TRAN_RECROSS;
                    log("Returning back to initial state, previous transition is a recrossing event\n");
                    recordFlag = true;
                }else{
                    jobStatus = ParallelReplicaJob::STATUS_NEWSTATE;
                    log("No recrossing, found new state\n");
                    *product = *current;
                    newStateStep = step; // remember the step when we are in a new state
                    if(parameters->parrepAutoStop){  // stop at transition; primarily for debugging
                        stopFlag = true;
                    }
                    recordFlag = false;
                }
            }
        }
         // new state was detected; determine refined transition time
        if(refineFlag && newStateFlag)
        {
            log("[Parallel Replica] Refining transition time\n");
            refFCalls = Potential::fcalls;
            refineStep = refine(mdBuffer, mdBufferLength, reactant);
            transitionStep = newStateStep - StateCheckInterval
                        - CorrSteps + refineStep*RecordInterval;
/* this is equivilant:
        transitionStep = transitionStep - StateCheckInterval
                         + refineStep*RecordInterval;
*/
            *transition = *mdBuffer[refineStep];
            refinedTime = (timeBuffer[refineStep]+timeBuffer[refineStep-1])/2.0;

            log("Found transition at step %ld; start decorrelation steps\n",
            transitionStep, parameters->parrepCorrTime);

            long corrBufferLength = long(CorrSteps/RecordInterval) + 1;
            if( (refineStep + corrBufferLength) < (mdBufferLength - 1) )
            {  // log("print here to debug %ld,%ld\n",refineStep+corrBufferLength,mdBufferLength);
                *product_relaxed = *mdBuffer[refineStep + corrBufferLength];
                corrTime = parameters->parrepCorrTime;
                transitionTime = refinedTime;
            }
            else
            {// log("nothing to say about this\n");
                *product_relaxed = *product;
                corrTime = time - refinedTime;
            }
            log("Correlation trajectory has been run for %.2f fs\n",corrTime*parameters->timeUnit);

            *transition_relaxed = *transition;
            relaxStatus = transition_relaxed->relax(true);
            if(!relaxStatus){
                jobStatus = ParallelReplicaJob::STATUS_BAD_RELAXFAILED;   
            }

            *product_relaxed = *product;
            relaxStatus = product_relaxed->relax(true);
            if(!relaxStatus){
                jobStatus = ParallelReplicaJob::STATUS_BAD_RELAXFAILED;
            }
      
            if(transition_relaxed->compare(product_relaxed)) {
                jobStatus = ParallelReplicaJob::STATUS_NEWSTATE;
            }else{
                jobStatus = ParallelReplicaJob::STATUS_NEWSTATE_CORR;
            }
            refineFCalls += Potential::fcalls - refFCalls;
            refineFlag = false;
        }

        // BOINC Progress
        if (step % 500 == 0) {
            // Since we only have a bundle size of 1 we can play with boinc_fraction_done
            // directly. When we have done parameters->mdSteps number of steps we aren't
            // quite done so I increase the max steps by 5%.
            boinc_fraction_done((double)step/(double)(parameters->mdSteps+0.05*parameters->mdSteps));
        }

        // stdout progress
        if ( (step % tenthSteps == 0) || (step == parameters->mdSteps) ) {
            double maxAtomDistance = current->perAtomNorm(*reactant);
            log("Progress: %3.0f%% ; Max displacement: %6.3lf ; Step: %7ld/%ld\n",
                (double)100.0*step/parameters->mdSteps, maxAtomDistance, step, parameters->mdSteps);
        }

        // we have run enough md steps; time to stop
        if (step == parameters->mdSteps-refineFCalls)
        {
            log("Achieved the specified MD simulation time\n");
            stopFlag = true;
        }else if (step > parameters->mdSteps-refineFCalls) {
            log("This trajectory will require more force calls\n");
            jobStatus = ParallelReplicaJob::STATUS_NEWSTATE_OVERFC;
            stopFlag = true;
        }
    }
    // calculate averages
    avgT = sumT/step;
    varT = sumT2/step - avgT*avgT;

    if (nboost > 0){
        log("\nTemperature : Average = %lf ; Stddev = %lf ; Factor = %lf; Boost = %lf\n\n",
        avgT, sqrt(varT), varT/avgT/avgT*nFreeCoord/2, sumboost/nboost);
    }else{
        log("\nTemperature : Average = %lf ; Stddev = %lf ; Factor = %lf\n\n",
        avgT, sqrt(varT), varT/avgT/avgT*nFreeCoord/2);
    }

    if (isfinite(avgT)==0)
    {
        log("Infinite average temperature, something went wrong!\n");
        newStateFlag = false;
        jobStatus = ParallelReplicaJob::STATUS_BAD_INFTEMP;
    }

    for(long i=0; i<mdBufferLength; i++){
        delete mdBuffer[i];
    }
    delete [] mdBuffer;
    delete [] timeBuffer;

    if(newStateFlag){
        return 1;
    }else{ 
        return 0;
    }
}

void ParallelReplicaJob::saveData(int state)
{
    FILE *fileResults, *fileReactant;

    std::string resultsFilename("results.dat");
    returnFiles.push_back(resultsFilename);

    fileResults = fopen(resultsFilename.c_str(), "wb");
    long totalFCalls = minimizeFCalls + mdFCalls + dephaseFCalls + refineFCalls+corrFCalls;

    fprintf(fileResults, "%s potential_type\n", parameters->potential.c_str());
    fprintf(fileResults, "%ld random_seed\n", parameters->randomSeed);
    fprintf(fileResults, "%lf potential_energy_reactant\n", reactant->getPotentialEnergy());
    fprintf(fileResults, "%ld total_force_calls\n", totalFCalls);
    fprintf(fileResults, "%ld force_calls_dephase\n", dephaseFCalls);
    fprintf(fileResults, "%ld force_calls_dynamics\n", mdFCalls);
    fprintf(fileResults, "%ld force_calls_minimize\n", minimizeFCalls);
    fprintf(fileResults, "%ld force_calls_refine\n", refineFCalls);
    fprintf(fileResults, "%ld force_calls_corr\n", corrFCalls);

    fprintf(fileResults, "%d termination_reason\n", jobStatus);
    fprintf(fileResults, "%d transition_found\n", (newStateFlag)?1:0);

    if(newStateFlag)
    {
        fprintf(fileResults, "%e transition_time_s\n", transitionTime*1.0e-15*parameters->timeUnit);
        fprintf(fileResults, "%e corr_time_s\n", corrTime*1.0e-15*parameters->timeUnit);
        fprintf(fileResults, "%lf potential_energy_product\n", product_relaxed->getPotentialEnergy());
        fprintf(fileResults, "%lf moved_distance\n",product_relaxed->distanceTo(*reactant));
    }

     
    fprintf(fileResults, "%e simulation_time_s\n", time*1.0e-15*parameters->timeUnit);
    fprintf(fileResults, "%lf speedup\n", time/(parameters->mdSteps*parameters->mdTimeStep));
    
    fclose(fileResults);

    std::string reactantFilename("reactant.con");
    returnFiles.push_back(reactantFilename);
    fileReactant = fopen(reactantFilename.c_str(), "wb");
    reactant->matter2con(fileReactant);
    fclose(fileReactant);

    if(newStateFlag)
    {
        FILE *fileProduct;
        std::string productFilename("product.con");
        returnFiles.push_back(productFilename);

        fileProduct = fopen(productFilename.c_str(), "wb");
        product_relaxed->matter2con(fileProduct);
        fclose(fileProduct);

        if(refineFCalls > 0)
        {
            FILE *fileTransition;
            std::string transitionFilename("transition.con");
            returnFiles.push_back(transitionFilename);

            fileTransition = fopen(transitionFilename.c_str(), "wb");
            transition->matter2con(fileTransition);
            fclose(fileTransition);
        }

        if(jobStatus == ParallelReplicaJob::STATUS_NEWSTATE_CORR)
        {
            FILE *fileMeta;
            std::string metaFilename("meta.con");
            returnFiles.push_back(metaFilename);
            fileMeta = fopen(metaFilename.c_str(), "wb");
            transition_relaxed->matter2con(fileMeta);
            fclose(fileMeta);
        }
    }
    return;
}


void ParallelReplicaJob::dephase()
{
    bool transitionFlag = false;
    long step, stepNew, loop;
    long dephaseSteps;
    long dephaseBufferLength, dephaseRefineStep;
    AtomMatrix velocity;
    Matter **dephaseBuffer;

    dephaseSteps = int(floor(parameters->parrepDephaseTime/parameters->mdTimeStep+0.5));
    Dynamics dephaseDynamics(current, parameters);
    log("Dephasing for %.2f fs\n",parameters->parrepDephaseTime*parameters->timeUnit);

    step = stepNew = loop = 0;

    //GH allocate dephase buffer memory
    dephaseBufferLength = dephaseSteps; // these are now the same variable
    dephaseBuffer = new Matter *[dephaseBufferLength];
    for(long i=0; i<dephaseBufferLength; i++){
        dephaseBuffer[i] = new Matter(parameters);
    }

    while(step < dephaseSteps)
    {
        loop++;
        // this should be allocated once, and of length dephaseSteps
        for(long i=step; i<dephaseBufferLength; i++)
        {
            dephaseDynamics.oneStep();
            *dephaseBuffer[i] = *current;
        }
    
        transitionFlag = checkState(current, reactant);

        if(transitionFlag)
        {
            dephaseRefineStep = refine(dephaseBuffer, dephaseBufferLength, reactant);
            log("Loop = %ld ; Dephase refine step = %ld\n", loop, dephaseRefineStep);
            transitionStep = dephaseRefineStep - 1; // check that this is correct
            transitionStep = (transitionStep > 0) ? transitionStep : 0;
            log("Dephasing warning: In a new state, invert the momentum and restart from step %ld\n", transitionStep);
            *current = *dephaseBuffer[transitionStep];
            velocity = current->getVelocities();
            velocity = velocity*(-1);
            current->setVelocities(velocity);
            step =  transitionStep;
        }
        else
        {
            step = dephaseBufferLength;
            log("Successful dephasing steps: %ld\n", step);
        }
        if( (parameters->parrepDephaseLoopStop) && (loop > parameters->parrepDephaseLoopMax) ) {
            log("Exceeded dephasing loop maximum; Dephasing steps: %ld\n", step);
            break;
        }
        log("Successful dephasing time: %.2f fs", step*parameters->mdTimeStep*parameters->timeUnit);

    }

    //GH deallocation dephase buffer memory
    for(long i=0; i<dephaseSteps; i++){
        delete dephaseBuffer[i];
    }
    delete [] dephaseBuffer;

}


bool ParallelReplicaJob::checkState(Matter *current, Matter *reactant)
{
    Matter tmp(parameters);
    tmp = *current;
    relaxStatus = tmp.relax(true);
    if(relaxStatus == false){
        jobStatus = ParallelReplicaJob::STATUS_BAD_RELAXFAILED;
    }
    if (tmp.compare(reactant)) {
        return false; 
    }
    return true;
}


long ParallelReplicaJob::refine(Matter *buff[], long length, Matter *reactant)
{
    //log("[Parallel Replica] Refining transition time.\n");

    bool midTest;
    long min, max, mid;

    min = 0;
    max = length - 1;

    while( (max-min) > 1 )
    {

        mid = min + (max-min)/2;
        midTest = checkState(buff[mid], reactant);

        if (midTest == false){
            min = mid;
        }
        else if (midTest == true){
            max = mid;
        }
        else {
            log("Refine step failed!\n");
            jobStatus = ParallelReplicaJob::STATUS_BAD_REFINEFAILED;
            exit(1);
        }
    }

    return (min+max)/2 + 1;
}
void ParallelReplicaJob::printEndStatus() {
    fprintf(stdout, "Final state: ");
    if(jobStatus == ParallelReplicaJob::STATUS_NEWSTATE)
        log("[ParallelReplica] New state found and is stable\n");

    else if(jobStatus == ParallelReplicaJob::STATUS_NEWSTATE_CORR) 
        log("[ParallelReplica] New state found; Correlated event detected; Metastable state has beeen saved as meta.con\n");  


    else if(jobStatus == ParallelReplicaJob::STATUS_NEWSTATE_OVERFC)
        log("[ParallelReplica] New state found; RefineFCalls exceeds the remaining MD steps after transition; Traj will take longer time to report\n");   
    else if(jobStatus == ParallelReplicaJob::STATUS_TRAN_NOTIME){
        log("[ParallelReplica] Insufficient force calls remaining to perform decorrelation and transition state refinement\n");
        log("[ParallelReplica] The last checkpoint that detected a transition will be reported\n");
    }

    else if(jobStatus == ParallelReplicaJob::STATUS_TRAN_RECROSS)
        log("[ParallelReplica] Transition found but has been determined as a recrossing event\n");

    else if(jobStatus == ParallelReplicaJob::STATUS_NOTRAN)
        log("[ParallelReplica] No transition was been found during the MD simulation time\n");

    else if( jobStatus == ParallelReplicaJob::STATUS_BAD_RELAXFAILED)
        log("[ParallelReplica] WARNING: Job failed in an optimization calculation\n");

    else if( jobStatus == ParallelReplicaJob::STATUS_BAD_REFINEFAILED)
        log("[ParallelReplica] WARNING: Job failed in the refinement process\n");

    else if( jobStatus == ParallelReplicaJob::STATUS_BAD_INFTEMP)
        log("[ParallelReplica] WARNING: Job running at INF temperature\n");

    else
        log("[ParallelReplica] Unknown jobStatus: %i!\n", jobStatus);
    return;
}
