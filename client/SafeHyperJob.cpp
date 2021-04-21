
#include <cstdlib>

#include "Matter.h"
#include "Dynamics.h"
#include "BondBoost.h"
#include "SafeHyperJob.h"
#include "Optimizer.h"
#include "Log.h"
#include <vector>

SafeHyperJob::SafeHyperJob(Parameters *params)
{
    parameters = params;
}

SafeHyperJob::~SafeHyperJob()
{
}

std::vector<std::string> SafeHyperJob::run(void)
{
    current = new Matter(parameters);
    reactant = new Matter(parameters);
    saddle = new Matter(parameters);
    product = new Matter(parameters);
    final = new Matter(parameters);
    final_tmp = new Matter(parameters);

    minimizeFCalls = mdFCalls = refineFCalls = dephaseFCalls = 0;
    time = 0.0;
    string reactantFilename = helper_functions::getRelevantFile(parameters->conFilename);
    current->con2matter(reactantFilename);

    log("\nMinimizing initial reactant\n");
    long refFCalls = Potential::fcalls;
    *reactant = *current;
    reactant->relax();
    minimizeFCalls += (Potential::fcalls - refFCalls);

    log("\nParallel Replica Dynamics, running\n\n");
    
    int status = dynamics();

    saveData(status);

    if(newStateFlag){
        log("Transition time: %.2e s\n", minCorrectedTime*1.0e-15*parameters->timeUnit);
    }else{
       log("No new state was found in %ld dynamics steps (%.3e s)\n",
           parameters->mdSteps, time*1.0e-15*parameters->timeUnit);
    }

    delete current;
    delete reactant;
    delete saddle;
    delete product;
    delete final_tmp;
    delete final;

    return returnFiles;
}

int SafeHyperJob::dynamics()
{
    bool transitionFlag = false, recordFlag = true, stopFlag = false, firstTransitFlag = false;
    long nFreeCoord = reactant->numberOfFreeAtoms()*3;
    long mdBufferLength, refFCalls;
    long step = 0, refineStep, newStateStep = 0; // check that newStateStep is set before used
    long nCheck = 0, nRecord = 0, nBoost = 0, nState = 0;
    long StateCheckInterval, RecordInterval;
    double kinE, kinT, avgT, varT;
    double kB = parameters->kB;
    double correctedTime = 0.0, sumCorrectedTime = 0.0, firstTransitionTime = 0.0;
    double Temp = 0.0, sumT = 0.0, sumT2 = 0.0; 
    double sumboost = 0.0, boost = 1.0, boostPotential = 0.0;
    double transitionTime_current=0.0, transitionTime_pre=0.0;
    AtomMatrix velocity;

    minCorrectedTime = 1.0e200;
    StateCheckInterval = int(parameters->parrepStateCheckInterval/parameters->mdTimeStep);
    RecordInterval = int(parameters->parrepRecordInterval/parameters->mdTimeStep);
    Temp = parameters->temperature;
    newStateFlag = metaStateFlag = false;

    mdBufferLength = long(StateCheckInterval/RecordInterval);
    Matter *mdBuffer[mdBufferLength];
    for(long i=0; i<mdBufferLength; i++) {
        mdBuffer[i] = new Matter(parameters);
    }
    timeBuffer = new double[mdBufferLength];
    biasBuffer = new double[mdBufferLength];

    Dynamics safeHyper(current, parameters);
    BondBoost bondBoost(current, parameters);

    if(parameters->biasPotential == Hyperdynamics::BOND_BOOST){
        bondBoost.initialize();
    }

    safeHyper.setThermalVelocity();

    // dephase the trajectory so that it is thermal and independent of others
    refFCalls = Potential::fcalls;
    dephase();
    dephaseFCalls = Potential::fcalls - refFCalls;

    log("\nStarting MD run\nTemperature: %.2f Kelvin\n"
        "Total Simulation Time: %.2f fs\nTime Step: %.2f fs\nTotal Steps: %ld\n\n", 
        Temp, 
        parameters->mdSteps*parameters->mdTimeStep*parameters->timeUnit,
        parameters->mdTimeStep*parameters->timeUnit,
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
            //log("step= %3d, boost = %10.5f",step,bondBoost.boost());
            boostPotential = bondBoost.boost();   
            boost = 1.0*exp(boostPotential/kB/Temp);
            time += parameters->mdTimeStep*boost;
            if (boost > 1.0){
                sumboost += boost;
                nBoost ++;
            }            
        } 

        kinE = current->getKineticEnergy();
        kinT = (2.0*kinE/nFreeCoord/kB); 
        sumT += kinT;
        sumT2 += kinT*kinT;
        //log("steps = %10d temp = %10.5f \n",step,kinT);

        safeHyper.oneStep();
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
                biasBuffer[nRecord] = boostPotential;
                nRecord++; // current location in the buffer
            }
        }

        // time to do a state check; if a transiton if found, stop recording
        if ( (nCheck == StateCheckInterval) && !newStateFlag )
        {
            nCheck = 0; // reinitialize check state counter
            nRecord = 0; // restart the buffer
            refFCalls = Potential::fcalls;
            transitionFlag = checkState(current, reactant);
            minimizeFCalls += Potential::fcalls - refFCalls;
            if(transitionFlag == true){
                nState ++;
                log("New State %ld: ",nState);
                *final_tmp = *current;
                transitionTime = time;
                newStateStep = step; // remember the step when we are in a new state
                transitionStep = newStateStep;
                firstTransitFlag = 1;
            }
        }
        //printf("step=%ld, time=%lf, biasPot=%lf\n",step,time,boostPotential);
        //Refine transition step
  
        if(transitionFlag)
        {
            //log("[Parallel Replica] Refining transition time.\n");
            refFCalls = Potential::fcalls;
            refineStep = refine(mdBuffer, mdBufferLength, reactant);

            transitionStep = newStateStep - StateCheckInterval + refineStep*RecordInterval;
            transitionTime_current = timeBuffer[refineStep];
            transitionTime = transitionTime_current - transitionTime_pre;
            transitionTime_pre = transitionTime_current;
            transitionPot = biasBuffer[refineStep];
            correctedTime = transitionTime * exp((-1)*transitionPot/kB/Temp);
            sumCorrectedTime += correctedTime;
            if ( nState == 1 ){
                firstTransitionTime = transitionTime;
            }
            //reverse the momenten;
            *current = *mdBuffer[refineStep-1];
            velocity = current->getVelocities();
            velocity = velocity*(-1);
            current->setVelocities(velocity);
            
            if (correctedTime < minCorrectedTime){
                minCorrectedTime = correctedTime;
                *saddle = *mdBuffer[refineStep];
                *final = *final_tmp;
            }
            log("tranisitonTime= %.3e s, biasPot= %.3f eV, correctedTime= %.3e s, sumCorrectedTime= %.3e s, minCorTime= %.3e s\n",transitionTime*1e-15*parameters->timeUnit,transitionPot,correctedTime*1e-15*parameters->timeUnit,sumCorrectedTime*1e-15*parameters->timeUnit, minCorrectedTime*1.0e-15*parameters->timeUnit);

            refineFCalls += Potential::fcalls - refFCalls;
            transitionFlag = false;
        }
    
            
        // we have run enough md steps; time to stop
        if (firstTransitFlag &&  sumCorrectedTime > firstTransitionTime) 
        {
            stopFlag = true;
            newStateFlag = true;
        }

        //stdout Progress
        if ( (step % tenthSteps == 0) || (step == parameters->mdSteps) ) {
            double maxAtomDistance = current->perAtomNorm(*reactant);
            log("progress: %3.0f%%, max displacement: %6.3lf, step %7ld/%ld\n",
                (double)100.0*step/parameters->mdSteps, maxAtomDistance, step, parameters->mdSteps);
        }
    }

    // calculate avearges
    avgT = sumT/step;
    varT = sumT2/step - avgT*avgT;
   
    if (nBoost > 0){
        log("\nTemperature : Average = %lf ; Stddev = %lf ; Factor = %lf; Boost = %lf\n\n",
        avgT, sqrt(varT), varT/avgT/avgT*nFreeCoord/2, sumboost/nBoost);
    }else{
        log("\nTemperature : Average = %lf ; Stddev = %lf ; Factor = %lf\n\n",
        avgT, sqrt(varT), varT/avgT/avgT*nFreeCoord/2);
    }
    if (isfinite(avgT)==0)
    {
        log("Infinite average temperature, something went wrong!\n");
        newStateFlag = false;
    }

    *product=*final; 
    // new state was detected; determine refined transition time
    for(long i=0; i<mdBufferLength; i++){
        delete mdBuffer[i];
    }
    delete [] timeBuffer;
    delete [] biasBuffer;

    if(newStateFlag){
        return 1;
    }else{
        return 0;
    }
}

void SafeHyperJob::saveData(int status)
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
        fprintf(fileResults, "%e transition_time_s\n", minCorrectedTime*1.0e-15*parameters->timeUnit);
        fprintf(fileResults, "%lf potential_energy_product\n", product->getPotentialEnergy());
        fprintf(fileResults, "%lf moved_distance\n",product->distanceTo(*reactant));
    }

     
    fprintf(fileResults, "%e simulation_time_s\n", time*1.0e-15*parameters->timeUnit);
    fprintf(fileResults, "%lf speedup\n", time/parameters->mdSteps/parameters->mdTimeStep);
    
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
    }
    return;
}


void SafeHyperJob::dephase()
{
    bool transitionFlag = false;
    long step, stepNew, loop;
    long DephaseSteps;
    long dephaseBufferLength, dephaseRefineStep;
    AtomMatrix velocity;

    DephaseSteps = int(parameters->parrepDephaseTime/parameters->mdTimeStep);
    Dynamics dephaseDynamics(current, parameters);
    log("Dephasing for %.2f fs\n",parameters->parrepDephaseTime*parameters->timeUnit);

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
        log("Successfully Dephased for %.2f fs", step*parameters->mdTimeStep*parameters->timeUnit);

    }
}


bool SafeHyperJob::checkState(Matter *current, Matter *reactant)
{
    Matter tmp(parameters);
    tmp = *current;
    tmp.relax(true);
    if (tmp.compare(reactant)) {
        return false; 
    }
    return true;
}


long SafeHyperJob::refine(Matter *buff[], long length, Matter *reactant)
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
            log("Refine step failed! \n");
            exit(1);
        }
    }

    return (min+max)/2 + 1;
}
