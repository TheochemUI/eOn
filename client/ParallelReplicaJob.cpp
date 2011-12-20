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
    final = new Matter(parameters);
    product = new Matter(parameters);

    minimizeFCalls = mdFCalls = refineFCalls = dephaseFCalls = 0;

    string reactant_passed = helper_functions::getRelevantFile(parameters->conFilename);
    current->con2matter(reactant_passed);

    log("\nMinimizing initial reactant\n");
    long refFCalls = Potential::fcalls;
    *reactant = *current;
    reactant->relax();
    minimizeFCalls += (Potential::fcalls - refFCalls);

    printf("Parallel Replica Dynamics, running\n\n");

    int status = dynamics();

    saveData(status);

    if(newStateFlag){
        printf("Transition time: %.2e s\n", transitionStep*1.018e-14);
    }else{
       printf("No new state was found in %ld dynamics steps (%.2f fs)\n",
           parameters->mdSteps, 10.18*parameters->mdSteps*parameters->mdTimeStep);
    }

    delete current;
    delete reactant;
    delete saddle;
    delete meta;
    delete final;
    delete product;

    return returnFiles;
}

int ParallelReplicaJob::dynamics()
{
    bool transitionFlag = false, recordFlag = true, stopFlag = false;
    long nFreeCoord = reactant->numberOfFreeAtoms()*3;
    long mdBufferLength, refFCalls;
    long step = 0, refineStep, productStep;
    long nCheck = 0, nRelax = 0, nRecord = 0;
    double time = 0.0;
    double kinE, kinT, avgT, varT,  kb = 1.0/11604.5;
    double sumT = 0.0, sumT2 = 0.0;

    newStateFlag = metaStateFlag = false;

    mdBufferLength = long(parameters->parrepStateCheckInterval/parameters->parrepRecordInterval) + 1;
    Matter *mdBuffer[mdBufferLength];
    for(long i=0; i<mdBufferLength; i++) {
        mdBuffer[i] = new Matter(parameters);
    }
    timeBuffer = new double[mdBufferLength];
    log("mdBufferLength = %ld\n", mdBufferLength);
    
    Dynamics parrepDynamics(current, parameters);
    BondBoost bondBoost(current, parameters);

    if(parameters->biasPotential == Hyperdynamics::BOND_BOOST){
        bondBoost.initial();
    }

    parrepDynamics.setThermalVelocity();
    dephase();

    printf("\nStarting MD run\nTemperature: %.2f Kelvin\n"
           "Total Time: %.2f fs\nTime Step: %.2f fs\n\n",
           parameters->temperature, 
           10.18*parameters->mdSteps*parameters->mdTimeStep,
           10.18*parameters->mdTimeStep);

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
            time += bondBoost.boost();
        } else {
            time += parameters->mdTimeStep;
        }

        kinE = current->getKineticEnergy();
        kinT = (2.*kinE/nFreeCoord/kb); 
        sumT += kinT;
        sumT2 += kinT*kinT;
       // printf("steps = %10d total_energy = %10.5f \n",nSteps,kinE+potE);

        parrepDynamics.oneStep();

        nCheck++; // count up to parameters->parrepStateCheckInterval before checking for a transition
        step++;

        // standard conditions; record mater object in the transition buffer
        if( parameters->parrepRefineTransition && recordFlag && !newStateFlag )
        {
            if( nCheck % parameters->parrepRecordInterval == 0 )
            {
                nRecord++; // current location in the buffer
                *mdBuffer[nRecord-1] = *current;
                timeBuffer[nRecord-1] = time;
            }
        }

//#ifndef NDEBUG
//        if (ncheck == parameters->parrepStateCheckInterval && !newState){
//           current->matter2xyz("movie", true);
//        }
//#endif

        // time to do a state check; if a transiton if found, stop recording
        if ( (nCheck == parameters->parrepStateCheckInterval) && !newStateFlag )
        {
            nCheck = 0; // reinitialize check state counter
            nRecord = 0; // restart the buffer
            transitionFlag = checkState(current, reactant);
            if(transitionFlag == true){
                recordFlag = false;
            }
        }

        // a transition has been detected; run an additional relaxSteps to see if the products are stable
        if ( transitionFlag && !newStateFlag )
        {
            nRelax++; // run relaxSteps before making a state check
            if (nRelax >= parameters->parrepRelaxSteps){
                // state check; reset counters
                nRelax = 0;
                nCheck = 0;
                nRecord = 0;
                newStateFlag = checkState(current, reactant);
                transitionFlag = false;
                if(newStateFlag == false){
                    recordFlag = true;
                }else{
                    printf("Found New State !\n");
                    *final = *current;
                    transitionStep = step; // remember the transition step
                    if(parameters->parrepAutoStop){  // stop at transition; primarily for debugging
                        stopFlag = true;
                    }
                    recordFlag = false;
                }
            }
        }

        // we have run enough md steps; time to stop
        if (step >= parameters->mdSteps){
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
            printf("progress: %3.0f%%, max displacement: %6.3lf, step %7ld/%ld\n",
                   (double)100.0*step/parameters->mdSteps, maxAtomDistance, step, parameters->mdSteps);
        }
    }
    avgT = sumT/step;
    varT = sumT2/step - avgT*avgT;

    log("\nTemperature : Average = %lf ; Variance = %lf ; Factor = %lf\n\n",
        avgT, varT, varT/avgT/avgT*nFreeCoord/2);

    if (isfinite(avgT)==0)
    {
        printf("Infinite average temperature, something went wrong!\n");
        newStateFlag = false;
    }

    productStep = transitionStep;

    // new state was detected; determine refined transition time
    if(parameters->parrepRefineTransition && newStateFlag)
    {
        refineStep = refine(mdBuffer, mdBufferLength, reactant);
        productStep = transitionStep - parameters->parrepStateCheckInterval -
                      parameters->parrepRelaxSteps + refineStep*parameters->parrepRecordInterval;
        *saddle = *mdBuffer[refineStep];
        transitionTime = timeBuffer[refineStep];

        printf("Found transition at step %ld, now running another %ld steps to allocate the product state\n",
            productStep, parameters->parrepRelaxSteps);

        long relaxBufferLength = int(parameters->parrepRelaxSteps/parameters->parrepRecordInterval) + 1;

        if(refineStep < (mdBufferLength - relaxBufferLength - 1) )
        {
            *final = *mdBuffer[productStep + relaxBufferLength];
        }
        else
        {
            // here, the final configuration should be obtained from the relax buffer -- fix this!
        }

        refFCalls = Potential::fcalls;
        *meta = *saddle;
        meta->relax(true);
        *product = *final;
        product->relax(true);
        minimizeFCalls += Potential::fcalls - refFCalls;

        if(*meta == *product){
           printf("Transition followed by a stable state.\n");
        }else{
           printf("Transition followed by a metastable state; product state taken after relaxation time.\n");
           transitionStep = transitionStep + parameters->parrepRelaxSteps;
           metaStateFlag = true;
        }
    }

    for(long i=0; i<mdBufferLength; i++){
        delete [] mdBuffer[i];
    }
    delete timeBuffer;

    if(newStateFlag){
        return 1;
    }else{
        return 0;
    }
}

void ParallelReplicaJob::saveData(int status)
{
    FILE *fileResults, *fileReactant, *fileProduct, *fileSaddle, *fileMeta;

    std::string resultsFilename("results.dat");
    returnFiles.push_back(resultsFilename);

    fileResults = fopen(resultsFilename.c_str(), "wb");
    long totalFCalls = minimizeFCalls + mdFCalls + dephaseFCalls + refineFCalls;

    fprintf(fileResults, "%d termination_reason\n", status);
    fprintf(fileResults, "%e transition_time_s\n", transitionStep*1.018e-14);
    fprintf(fileResults, "%ld random_seed\n", parameters->randomSeed);
    fprintf(fileResults, "%lf potential_energy_reactant\n", reactant->getPotentialEnergy());
    fprintf(fileResults, "%lf potential_energy_product\n", product->getPotentialEnergy());
    fprintf(fileResults, "%s potential_type\n", parameters->potential.c_str());
    fprintf(fileResults, "%ld total_force_calls\n", totalFCalls);
    fprintf(fileResults, "%ld force_calls_dephase\n", dephaseFCalls);
    fprintf(fileResults, "%ld force_calls_dynamics\n", mdFCalls);
    fprintf(fileResults, "%ld force_calls_minimize\n", minimizeFCalls);
    fprintf(fileResults, "%ld force_calls_refine\n", refineFCalls);

    fprintf(fileResults, "%lf moved_distance\n",final->distanceTo(*reactant));
    fclose(fileResults);

    std::string reactantFilename("reactant.con");
    returnFiles.push_back(reactantFilename);
    fileReactant = fopen(reactantFilename.c_str(), "wb");
    reactant->matter2con(fileReactant);
    fclose(fileReactant);

    std::string productFilename("product.con");
    returnFiles.push_back(productFilename);

    fileProduct = fopen(productFilename.c_str(), "wb");
    product->matter2con(fileProduct);
    fclose(fileProduct);

    std::string saddleFilename("saddle.con");
    returnFiles.push_back(saddleFilename);

    fileSaddle = fopen(saddleFilename.c_str(), "wb");
    saddle->matter2con(fileSaddle);
    fclose(fileSaddle);

    if(metaStateFlag){
        std::string metaFilename("meta.con");
        returnFiles.push_back(metaFilename);
        fileMeta = fopen(metaFilename.c_str(), "wb");
        meta->matter2con(fileMeta);
        fclose(fileMeta);
    }
    return;
}


void ParallelReplicaJob::dephase()
{
    long step, stepNew, loop;
    long dephaseBufferLength, dephaseRefineStep;
    long refFCalls;
    AtomMatrix velocity;

    Dynamics dephaseDynamics(current, parameters);
    printf("Dephasing for %ld steps\n", parameters->parrepDephaseSteps);

    step = stepNew = loop = 0;

    while(step < parameters->parrepDephaseSteps)
    {
        // this should be allocated once, and of length parameters->parrepDephaseSteps
        dephaseBufferLength = parameters->parrepDephaseSteps - step;
        loop++;
        Matter *dephaseBuffer[dephaseBufferLength];

        refFCalls = Potential::fcalls;
        for(long i=0; i<dephaseBufferLength; i++)
        {
           dephaseBuffer[i] = new Matter(parameters);
           dephaseDynamics.oneStep();
           *dephaseBuffer[i] = *current;
        }
        dephaseFCalls += Potential::fcalls - refFCalls;
    
        newStateFlag = checkState(current, reactant);

        if(newStateFlag)
        {
            dephaseRefineStep = refine(dephaseBuffer, dephaseBufferLength, reactant);
            printf("loop = %ld; dephase refine step = %ld\n", loop, dephaseRefineStep);
            transitionStep = dephaseRefineStep - 1; // check that this is correct
            transitionStep = (transitionStep > 0) ? transitionStep : 0;
            printf("Dephasing warning: in a new state, inverse the momentum and restart from step %ld\n", step+transitionStep);
            *current = *dephaseBuffer[transitionStep];
            velocity = current->getVelocities();
            velocity = velocity*(-1);
            current->setVelocities(velocity);
            step = step + transitionStep;
        } else {
            step = step + dephaseBufferLength;
            printf("Successful dephasing for %ld steps \n", step);
        }

        for(long i=0; i<dephaseBufferLength; i++)
        {
           delete dephaseBuffer[i];
        }

        if( (parameters->parrepDephaseLoopStop) && (loop > parameters->parrepDephaseLoopMax) ) {
            printf("Reach dephase loop maximum, stop dephasing! Dephased for %ld steps\n ", step);
            break;
        }
    }
}


bool ParallelReplicaJob::checkState(Matter *current, Matter *reactant)
{
    Matter tmp(parameters);
    tmp = *current;

    long refFCalls = Potential::fcalls;
    tmp.relax(true);
    minimizeFCalls += Potential::fcalls - refFCalls;

    if (tmp == *reactant) {
        return false;
    }
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
