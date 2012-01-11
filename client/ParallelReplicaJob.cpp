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

ParallelReplicaJob::ParallelReplicaJob(Parameters *parameters_passed)
{
    parameters = parameters_passed;
}

ParallelReplicaJob::~ParallelReplicaJob()
{
}

std::vector<std::string> ParallelReplicaJob::run(void)
{
    current = new Matter(parameters);
    reactant = new Matter(parameters);
    saddle = new Matter(parameters);
    meta = new Matter(parameters);
    product = new Matter(parameters);
    final = new Matter(parameters);

    minimizeFCalls = mdFCalls = refineFCalls = dephaseFCalls = 0;
    time = 0.0;
    string reactant_passed = helper_functions::getRelevantFile(parameters->conFilename);
    current->con2matter(reactant_passed);

    log("\nMinimizing initial reactant\n");
    long refFCalls = Potential::fcalls;
    *reactant = *current;
    reactant->relax();
    minimizeFCalls += (Potential::fcalls - refFCalls);

    log("\nParallel Replica Dynamics, running\n\n");
    
    int status = dynamics();

    saveData(status);

    if(newStateFlag){
        log("Transition time: %.2e s\n", transitionTime);
    }else{
       log("No new state was found in %ld dynamics steps (%.2f fs)\n",
           parameters->mdSteps, time);
    }

    delete current;
    delete reactant;
    delete saddle;
    delete meta;
    delete product;
    delete final;

    return returnFiles;
}

int ParallelReplicaJob::dynamics()
{
    bool transitionFlag = false, recordFlag = true, stopFlag = false;
    long nFreeCoord = reactant->numberOfFreeAtoms()*3;
    long mdBufferLength, refFCalls;
    long step = 0, refineStep, newStateStep = 0; // check that newStateStep is set before used
    long nCheck = 0, nRelax = 0, nRecord = 0;
    long StateCheckInterval, RecordInterval, RelaxSteps;
    double kinE, kinT, avgT, varT,  kb = 1.0/11604.5;
    double sumT = 0.0, sumT2 = 0.0;

    StateCheckInterval = int(parameters->parrepStateCheckInterval/parameters->mdTimeStepInput);
    RecordInterval = int(parameters->parrepRecordInterval/parameters->mdTimeStepInput);
    RelaxSteps = int(parameters->parrepRelaxTime/parameters->mdTimeStepInput);
    newStateFlag = metaStateFlag = false;

    mdBufferLength = long(StateCheckInterval/RecordInterval);
    Matter *mdBuffer[mdBufferLength];
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

    log("\nStarting MD run\nTemperature: %.2f Kelvin\n"
        "Total Simulation Time: %.2f fs\nTime Step: %.2f fs\nTotal Steps: %ld\n\n", 
        parameters->temperature, 
        parameters->mdSteps*parameters->mdTimeStepInput,
        parameters->mdTimeStepInput,
        parameters->mdSteps);
    log("MD buffer length: %ld\n", mdBufferLength);

    long tenthSteps = parameters->mdSteps/10;
    //This prevents and edge case division by zero if mdSteps is < 10
    if (tenthSteps == 0) {
        tenthSteps = parameters->mdSteps;
    }

    // loop dynamics iterations until some condition tells us to stop
    while(!stopFlag)
    {
        if( (parameters->biasPotential == Hyperdynamics::BOND_BOOST) && !newStateFlag ) {
            // GH: boost should be a unitless factor, multipled by TimeStep to get the boosted time
            time += parameters->mdTimeStepInput*bondBoost.boost();
        } else {
            time += parameters->mdTimeStepInput;
        }

        kinE = current->getKineticEnergy();
        kinT = (2.*kinE/nFreeCoord/kb); 
        sumT += kinT;
        sumT2 += kinT*kinT;
       // log("steps = %10d total_energy = %10.5f \n",nSteps,kinE+potE);

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
        if ( (nCheck == StateCheckInterval) && !newStateFlag )
        {
            nCheck = 0; // reinitialize check state counter
            nRecord = 0; // restart the buffer
            refFCalls = Potential::fcalls;
            transitionFlag = checkState(current, reactant);
            minimizeFCalls += Potential::fcalls - refFCalls;
            if(transitionFlag == true){
                transitionStep = step;
                recordFlag = false;
            }
        }

        // a transition has been detected; run an additional relaxSteps to see if the products are stable
        if ( transitionFlag && !newStateFlag )
        {
            nRelax++; // run relaxSteps before making a state check
            if (nRelax > RelaxSteps){
                // state check; reset counters
                nRelax = 0;
                nCheck = 0;
                nRecord = 0;
                refFCalls = Potential::fcalls;
                newStateFlag = checkState(current, reactant);
                minimizeFCalls += Potential::fcalls - refFCalls;
                transitionFlag = false;
                if(newStateFlag == false){
                    recordFlag = true;
                }else{
                    log("Found New State.\n");
                    *final = *current;
                    //refFCalls = Potential::fcalls;
                    //product->relax(true);
                   // minimizeFCalls += Potential::fcalls - refFCalls;
                    newStateStep = step; // remember the step when we are in a new state
                    if(parameters->parrepAutoStop){  // stop at transition; primarily for debugging
                        stopFlag = true;
                    }
                    recordFlag = false;
                }
            }
        }

        // we have run enough md steps; time to stop
        if (step >= parameters->mdSteps)
        {
            stopFlag = true;
        }

        //BOINC Progress
        if (step % 500 == 0) {
            // Since we only have a bundle size of 1 we can play with boinc_fraction_done
            // directly. When we have done parameters->mdSteps number of steps we aren't
            // quite done so I increase the max steps by 5%.
            boinc_fraction_done((double)step/(double)(parameters->mdSteps+0.05*parameters->mdSteps));
        }

        //stdout Progress
        if ( (step % tenthSteps == 0) || (step == parameters->mdSteps) ) {
            double maxAtomDistance = current->perAtomNorm(*reactant);
            log("progress: %3.0f%%, max displacement: %6.3lf, step %7ld/%ld\n",
                (double)100.0*step/parameters->mdSteps, maxAtomDistance, step, parameters->mdSteps);
        }
    }

    // calculate avearges
    avgT = sumT/parameters->mdSteps;
    varT = sumT2/parameters->mdSteps - avgT*avgT;

    log("\nTemperature : Average = %lf ; Stddev = %lf ; Factor = %lf\n\n",
        avgT, sqrt(varT), varT/avgT/avgT*nFreeCoord/2);

    if (isfinite(avgT)==0)
    {
        log("Infinite average temperature, something went wrong!\n");
        newStateFlag = false;
    }

    // new state was detected; determine refined transition time
    if(parameters->parrepRefineTransition && newStateFlag)
    {
        refFCalls = Potential::fcalls;
        refineStep = refine(mdBuffer, mdBufferLength, reactant);

        transitionStep = newStateStep - StateCheckInterval
                        - RelaxSteps + refineStep*RecordInterval;
/* this is equivilant:
        transitionStep = transitionStep - StateCheckInterval
                         + refineStep*RecordInterval;
*/
        *saddle = *mdBuffer[refineStep];
        transitionTime = timeBuffer[refineStep];

        log("Found transition at step %ld, now running another %.2f fs to allocate the product state.\n",
            transitionStep, parameters->parrepRelaxTime);

        long relaxBufferLength = long(RelaxSteps/RecordInterval) + 1;

        if( (refineStep + relaxBufferLength) < (mdBufferLength - 1) )
        {
           // log("print here to debug %ld,%ld\n",refineStep+relaxBufferLength,mdBufferLength);
            *product = *mdBuffer[refineStep + relaxBufferLength];
        }
        else
        {
           // log("nothing to say about this\n");
            *product = *final;
            // here, the final configuration should be obtained from the relax buffer -- fix this!
        }

        *meta = *saddle;
        meta->relax(true);
        product->relax(true);

        if(*meta == *product){
           log("Transition followed by a stable state.\n");
        }else{
           log("Transition followed by a metastable state; product state taken after relaxation time.\n");
           transitionStep = transitionStep + RelaxSteps;
           transitionTime = transitionTime + parameters->parrepRelaxTime;
           metaStateFlag = true;
        }
        refineFCalls += Potential::fcalls - refFCalls;
    }

    for(long i=0; i<mdBufferLength; i++){
        delete mdBuffer[i];
    }
    delete [] timeBuffer;

    if(newStateFlag){
        return 1;
    }else{
        return 0;
    }
}

void ParallelReplicaJob::saveData(int status)
{
    FILE *fileResults, *fileReactant;

    std::string resultsFilename("results.dat");
    returnFiles.push_back(resultsFilename);

    fileResults = fopen(resultsFilename.c_str(), "wb");
    long totalFCalls = minimizeFCalls + mdFCalls + dephaseFCalls + refineFCalls;

    fprintf(fileResults, "%s potential_type\n", parameters->potential.c_str());
    fprintf(fileResults, "%ld random_seed\n", parameters->randomSeed);
    fprintf(fileResults, "%lf potential_energy_reactant\n", reactant->getPotentialEnergy());
    fprintf(fileResults, "%ld total_force_calls\n", totalFCalls);
    fprintf(fileResults, "%ld force_calls_dephase\n", dephaseFCalls);
    fprintf(fileResults, "%ld force_calls_dynamics\n", mdFCalls);
    fprintf(fileResults, "%ld force_calls_minimize\n", minimizeFCalls);
    fprintf(fileResults, "%ld force_calls_refine\n", refineFCalls);
 

//    fprintf(fileResults, "%d termination_reason\n", status);
    fprintf(fileResults, "%d transition_found\n", (newStateFlag)?1:0);

    if(newStateFlag)
    {
        fprintf(fileResults, "%e transition_time_s\n", transitionTime);
        fprintf(fileResults, "%lf potential_energy_product\n", product->getPotentialEnergy());
        fprintf(fileResults, "%lf moved_distance\n",product->distanceTo(*reactant));
    }
    else
    { 
        fprintf(fileResults, "%e simulation_time_s\n", time);
    }
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
        product->matter2con(fileProduct);
        fclose(fileProduct);

        if(parameters->parrepRefineTransition)
        {
            FILE *fileSaddle;
            std::string saddleFilename("saddle.con");
            returnFiles.push_back(saddleFilename);

            fileSaddle = fopen(saddleFilename.c_str(), "wb");
            saddle->matter2con(fileSaddle);
            fclose(fileSaddle);
        }

        if(metaStateFlag)
        {
            FILE *fileMeta;
            std::string metaFilename("meta.con");
            returnFiles.push_back(metaFilename);
            fileMeta = fopen(metaFilename.c_str(), "wb");
            meta->matter2con(fileMeta);
            fclose(fileMeta);
        }
    }
    return;
}


void ParallelReplicaJob::dephase()
{
    bool transitionFlag = false;
    long step, stepNew, loop;
    long DephaseSteps;
    long dephaseBufferLength, dephaseRefineStep;
    AtomMatrix velocity;

    DephaseSteps = int(parameters->parrepDephaseTime/parameters->mdTimeStepInput);
    Dynamics dephaseDynamics(current, parameters);
    log("Dephasing for %.2f fs\n",parameters->parrepDephaseTime);

    step = stepNew = loop = 0;

    while(step < DephaseSteps)
    {
        // this should be allocated once, and of length DephaseSteps
        dephaseBufferLength = DephaseSteps - step;
        loop++;
        Matter *dephaseBuffer[dephaseBufferLength];

        for(long i=0; i<dephaseBufferLength; i++)
        {
            dephaseBuffer[i] = new Matter(parameters);
            dephaseDynamics.oneStep();
            *dephaseBuffer[i] = *current;
        }
    
        transitionFlag = checkState(current, reactant);

        if(transitionFlag)
        {
            dephaseRefineStep = refine(dephaseBuffer, dephaseBufferLength, reactant);
            log("loop = %ld; dephase refine step = %ld\n", loop, dephaseRefineStep);
            transitionStep = dephaseRefineStep - 1; // check that this is correct
            transitionStep = (transitionStep > 0) ? transitionStep : 0;
            log("Dephasing warning: in a new state, inverse the momentum and restart from step %ld\n", step+transitionStep);
            *current = *dephaseBuffer[transitionStep];
            velocity = current->getVelocities();
            velocity = velocity*(-1);
            current->setVelocities(velocity);
            step = step + transitionStep;
        }
        else
        {
            step = step + dephaseBufferLength;
            //log("Successful dephasing for %.2f steps \n", step);
        }

        for(long i=0; i<dephaseBufferLength; i++)
        {
           delete dephaseBuffer[i];
        }

        if( (parameters->parrepDephaseLoopStop) && (loop > parameters->parrepDephaseLoopMax) ) {
            log("Reach dephase loop maximum, stop dephasing! Dephased for %ld steps\n ", step);
            break;
        }
        log("Successfully Dephased for %.2f fs", step*parameters->mdTimeStepInput);

    }
}


bool ParallelReplicaJob::checkState(Matter *current, Matter *reactant)
{
    Matter tmp(parameters);
    tmp = *current;
    tmp.relax(true);
    if (tmp == *reactant) {
        return false; }
    return true;
}


long ParallelReplicaJob::refine(Matter *buff[], long length, Matter *reactant)
{
    log("[Parallel Replica] Refining transition time.\n");

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
            log("Refine step failed! \n");
            exit(1);
        }
    }

    return (min+max)/2 + 1;
}
