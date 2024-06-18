#include "SafeHyperJob.h"
#include "BondBoost.h"
#include "Dynamics.h"
#include "Optimizer.h"

std::vector<std::string> SafeHyperJob::run(void) {
  // TODO: Rework
  current = new Matter(pot, params);
  reactant = new Matter(pot, params);
  saddle = new Matter(pot, params);
  product = new Matter(pot, params);
  final_img = new Matter(pot, params);
  final_img_tmp = new Matter(pot, params);

  minimizeFCalls = mdFCalls = refineFCalls = dephaseFCalls = 0;
  time = 0.0;
  string reactantFilename =
      helper_functions::getRelevantFile(params->conFilename);
  current->con2matter(reactantFilename);

  SPDLOG_LOGGER_DEBUG(log, "Minimizing initial reactant");
  long refFCalls = Potential::fcalls;
  *reactant = *current;
  reactant->relax();
  // minimizeFCalls += (Potential::fcalls - refFCalls);

  SPDLOG_LOGGER_DEBUG(log, "Parallel Replica Dynamics, running");

  int status = dynamics();

  saveData(status);

  if (newStateFlag) {
    SPDLOG_LOGGER_DEBUG(log, "Transition time: {:.2e} s",
                        minCorrectedTime * 1.0e-15 * params->timeUnit);
  } else {
    SPDLOG_LOGGER_DEBUG(
        log, "No new state was found in {} dynamics steps ({:.3e} s)",
        params->mdSteps, time * 1.0e-15 * params->timeUnit);
  }

  delete current;
  delete reactant;
  delete saddle;
  delete product;
  delete final_img_tmp;
  delete final_img;

  return returnFiles;
}

int SafeHyperJob::dynamics() {
  bool transitionFlag = false, recordFlag = true, stopFlag = false,
       firstTransitFlag = false;
  long nFreeCoord = reactant->numberOfFreeAtoms() * 3;
  long mdBufferLength, refFCalls;
  long step = 0, refineStep,
       newStateStep = 0; // check that newStateStep is set before used
  long nCheck = 0, nRecord = 0, nBoost = 0, nState = 0;
  long StateCheckInterval, RecordInterval;
  double kinE, kinT, avgT, varT;
  double kB = params->kB;
  double correctedTime = 0.0, sumCorrectedTime = 0.0, firstTransitionTime = 0.0;
  double Temp = 0.0, sumT = 0.0, sumT2 = 0.0;
  double sumboost = 0.0, boost = 1.0, boostPotential = 0.0;
  double transitionTime_current = 0.0, transitionTime_pre = 0.0;
  AtomMatrix velocity;

  minCorrectedTime = 1.0e200;
  StateCheckInterval =
      int(params->parrepStateCheckInterval / params->mdTimeStep);
  RecordInterval = int(params->parrepRecordInterval / params->mdTimeStep);
  Temp = params->temperature;
  newStateFlag = metaStateFlag = false;

  mdBufferLength = long(StateCheckInterval / RecordInterval);
  Matter *mdBuffer[mdBufferLength];
  for (long i = 0; i < mdBufferLength; i++) {
    mdBuffer[i] = new Matter(pot, params);
  }
  timeBuffer = new double[mdBufferLength];
  biasBuffer = new double[mdBufferLength];

  Dynamics safeHyper(current, params.get());
  BondBoost bondBoost(current, params.get());

  if (params->biasPotential == Hyperdynamics::BOND_BOOST) {
    bondBoost.initialize();
  }

  safeHyper.setThermalVelocity();

  // dephase the trajectory so that it is thermal and independent of others
  // refFCalls = Potential::fcalls;
  dephase();
  // dephaseFCalls = Potential::fcalls - refFCalls;

  SPDLOG_LOGGER_DEBUG(
      log,
      "Starting MD run\nTemperature: {:.2f} Kelvin\n"
      "Total Simulation Time: {:.2f} fs\nTime Step: {:.2f} fs\nTotal Steps: {}",
      Temp, params->mdSteps * params->mdTimeStep * params->timeUnit,
      params->mdTimeStep * params->timeUnit, params->mdSteps);
  SPDLOG_LOGGER_DEBUG(log, "MD buffer length: {}", mdBufferLength);

  long tenthSteps = params->mdSteps / 10;
  // This prevents and edge case division by zero if mdSteps is < 10
  if (tenthSteps == 0) {
    tenthSteps = params->mdSteps;
  }

  // loop dynamics iterations until some condition tells us to stop
  while (!stopFlag) {
    if ((params->biasPotential == Hyperdynamics::BOND_BOOST) && !newStateFlag) {
      // GH: boost should be a unitless factor, multipled by TimeStep to get the
      // boosted time
      // SPDLOG_LOGGER_DEBUG(log, "step= {} , boost = {:.5f}", step,
      // bondBoost.boost());
      boostPotential = bondBoost.boost();
      boost = 1.0 * exp(boostPotential / kB / Temp);
      time += params->mdTimeStep * boost;
      if (boost > 1.0) {
        sumboost += boost;
        nBoost++;
      }
    }

    kinE = current->getKineticEnergy();
    kinT = (2.0 * kinE / nFreeCoord / kB);
    sumT += kinT;
    sumT2 += kinT * kinT;
    // SPDLOG_LOGGER_DEBUG(log, "steps = {:10} temp = {:10.5f}", step, kinT);

    safeHyper.oneStep();
    mdFCalls++;

    nCheck++; // count up to params->parrepStateCheckInterval before checking
              // for a transition
    step++;
    // SPDLOG_LOGGER_DEBUG(log, "step = {:4}, time = {:10.4f}", step, time);
    //  standard conditions; record mater object in the transition buffer
    if (params->parrepRefineTransition && recordFlag && !newStateFlag) {
      if (nCheck % RecordInterval == 0) {
        *mdBuffer[nRecord] = *current;
        timeBuffer[nRecord] = time;
        biasBuffer[nRecord] = boostPotential;
        nRecord++; // current location in the buffer
      }
    }

    // time to do a state check; if a transiton if found, stop recording
    if ((nCheck == StateCheckInterval) && !newStateFlag) {
      nCheck = 0;  // reinitialize check state counter
      nRecord = 0; // restart the buffer
      refFCalls = Potential::fcalls;
      transitionFlag = checkState(current, reactant);
      // minimizeFCalls += Potential::fcalls - refFCalls;
      if (transitionFlag == true) {
        nState++;
        SPDLOG_LOGGER_DEBUG(log, "New State {}: ", nState);
        *final_img_tmp = *current;
        transitionTime = time;
        newStateStep = step; // remember the step when we are in a new state
        transitionStep = newStateStep;
        firstTransitFlag = 1;
      }
    }
    // printf("step=%ld, time=%lf, biasPot=%lf\n",step,time,boostPotential);
    // Refine transition step

    if (transitionFlag) {
      // SPDLOG_LOGGER_DEBUG(log, "[Parallel Replica] Refining transition
      // time.");
      refFCalls = Potential::fcalls;
      refineStep = refine(mdBuffer, mdBufferLength, reactant);

      transitionStep =
          newStateStep - StateCheckInterval + refineStep * RecordInterval;
      transitionTime_current = timeBuffer[refineStep];
      transitionTime = transitionTime_current - transitionTime_pre;
      transitionTime_pre = transitionTime_current;
      transitionPot = biasBuffer[refineStep];
      correctedTime = transitionTime * exp((-1) * transitionPot / kB / Temp);
      sumCorrectedTime += correctedTime;
      if (nState == 1) {
        firstTransitionTime = transitionTime;
      }
      // reverse the momenten;
      *current = *mdBuffer[refineStep - 1];
      velocity = current->getVelocities();
      velocity = velocity * (-1);
      current->setVelocities(velocity);

      if (correctedTime < minCorrectedTime) {
        minCorrectedTime = correctedTime;
        *saddle = *mdBuffer[refineStep];
        *final_img = *final_img_tmp;
      }
      SPDLOG_LOGGER_DEBUG(log,
                          "tranisitonTime= {:.3e} s, biasPot= {:.3f} eV, "
                          "correctedTime= {:.3e} s, "
                          "sumCorrectedTime= {:.3e} s, minCorTime= {:.3e} s",
                          transitionTime * 1e-15 * params->timeUnit,
                          transitionPot,
                          correctedTime * 1e-15 * params->timeUnit,
                          sumCorrectedTime * 1e-15 * params->timeUnit,
                          minCorrectedTime * 1.0e-15 * params->timeUnit);

      // refineFCalls += Potential::fcalls - refFCalls;
      transitionFlag = false;
    }

    // we have run enough md steps; time to stop
    if (firstTransitFlag && sumCorrectedTime > firstTransitionTime) {
      stopFlag = true;
      newStateFlag = true;
    }

    // stdout Progress
    if ((step % tenthSteps == 0) || (step == params->mdSteps)) {
      double maxAtomDistance = current->perAtomNorm(*reactant);
      SPDLOG_LOGGER_DEBUG(
          log, "progress: {:.0f}%, max displacement: {:.3f}, step {} / {}",
          static_cast<double>(100.0 * step / params->mdSteps), maxAtomDistance,
          step, params->mdSteps);
    }
  }

  // calculate averages
  avgT = sumT / step;
  varT = sumT2 / step - avgT * avgT;

  if (nBoost > 0) {
    SPDLOG_LOGGER_DEBUG(log,
                        "Temperature : Average = {:.6f} ; Stddev = {:.6f} ; "
                        "Factor = {:.6f}; Boost = {:.6f}",
                        avgT, sqrt(varT), varT / avgT / avgT * nFreeCoord / 2,
                        sumboost / nBoost);
  } else {
    SPDLOG_LOGGER_DEBUG(
        log,
        "Temperature : Average = {:.6f} ; Stddev = {:.6f} ; Factor = {:.6f}",
        avgT, sqrt(varT), varT / avgT / avgT * nFreeCoord / 2);
  }
  if (std::isfinite(avgT) == 0) {
    SPDLOG_LOGGER_DEBUG(log,
                        "Infinite average temperature, something went wrong!");
    newStateFlag = false;
  }

  *product = *final_img;
  // new state was detected; determine refined transition time
  for (long i = 0; i < mdBufferLength; i++) {
    delete mdBuffer[i];
  }
  delete[] timeBuffer;
  delete[] biasBuffer;

  if (newStateFlag) {
    return 1;
  } else {
    return 0;
  }
}

void SafeHyperJob::saveData(int status) {
  FILE *fileResults, *fileReactant;

  std::string resultsFilename("results.dat");
  returnFiles.push_back(resultsFilename);

  fileResults = fopen(resultsFilename.c_str(), "wb");
  // long totalFCalls = minimizeFCalls + mdFCalls + dephaseFCalls +
  // refineFCalls;

  fprintf(
      fileResults, "%s potential_type\n",
      std::string{magic_enum::enum_name<PotType>(params->potential)}.c_str());
  fprintf(fileResults, "%ld random_seed\n", params->randomSeed);
  fprintf(fileResults, "%lf potential_energy_reactant\n",
          reactant->getPotentialEnergy());
  // fprintf(fileResults, "%ld total_force_calls\n", totalFCalls);
  fprintf(fileResults, "%ld force_calls_dephase\n", dephaseFCalls);
  fprintf(fileResults, "%ld force_calls_dynamics\n", mdFCalls);
  fprintf(fileResults, "%ld force_calls_minimize\n", minimizeFCalls);
  // fprintf(fileResults, "%ld force_calls_refine\n", refineFCalls);

  //    fprintf(fileResults, "%d termination_reason\n", status);
  fprintf(fileResults, "%d transition_found\n", (newStateFlag) ? 1 : 0);

  if (newStateFlag) {
    fprintf(fileResults, "%e transition_time_s\n",
            minCorrectedTime * 1.0e-15 * params->timeUnit);
    fprintf(fileResults, "%lf potential_energy_product\n",
            product->getPotentialEnergy());
    fprintf(fileResults, "%lf moved_distance\n",
            product->distanceTo(*reactant));
  }

  fprintf(fileResults, "%e simulation_time_s\n",
          time * 1.0e-15 * params->timeUnit);
  fprintf(fileResults, "%lf speedup\n",
          time / params->mdSteps / params->mdTimeStep);

  fclose(fileResults);

  std::string reactantFilename("reactant.con");
  returnFiles.push_back(reactantFilename);
  fileReactant = fopen(reactantFilename.c_str(), "wb");
  reactant->matter2con(fileReactant);
  fclose(fileReactant);

  if (newStateFlag) {
    FILE *fileProduct;
    std::string productFilename("product.con");
    returnFiles.push_back(productFilename);

    fileProduct = fopen(productFilename.c_str(), "wb");
    product->matter2con(fileProduct);
    fclose(fileProduct);

    if (params->parrepRefineTransition) {
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

void SafeHyperJob::dephase() {
  bool transitionFlag = false;
  long step, stepNew, loop;
  long DephaseSteps;
  long dephaseBufferLength, dephaseRefineStep;
  AtomMatrix velocity;

  DephaseSteps = int(params->parrepDephaseTime / params->mdTimeStep);
  Dynamics dephaseDynamics(current, params.get());
  SPDLOG_LOGGER_DEBUG(log, "Dephasing for {:.2f} fs",
                      params->parrepDephaseTime * params->timeUnit);

  step = stepNew = loop = 0;

  while (step < DephaseSteps) {
    // this should be allocated once, and of length DephaseSteps
    dephaseBufferLength = DephaseSteps - step;
    loop++;
    Matter *dephaseBuffer[dephaseBufferLength];

    for (long i = 0; i < dephaseBufferLength; i++) {
      dephaseBuffer[i] = new Matter(pot, params);
      dephaseDynamics.oneStep();
      *dephaseBuffer[i] = *current;
    }

    transitionFlag = checkState(current, reactant);

    if (transitionFlag) {
      dephaseRefineStep = refine(dephaseBuffer, dephaseBufferLength, reactant);
      SPDLOG_LOGGER_DEBUG(log, "loop = {}; dephase refine step = {}", loop,
                          dephaseRefineStep);
      transitionStep = dephaseRefineStep - 1; // check that this is correct
      transitionStep = (transitionStep > 0) ? transitionStep : 0;
      SPDLOG_LOGGER_DEBUG(log,
                          "Dephasing warning: in a new state, inverse the "
                          "momentum and restart "
                          "from step {}",
                          step + transitionStep);
      *current = *dephaseBuffer[transitionStep];
      velocity = current->getVelocities();
      velocity = velocity * (-1);
      current->setVelocities(velocity);
      step = step + transitionStep;
    } else {
      step = step + dephaseBufferLength;
      // SPDLOG_LOGGER_DEBUG(log, "Successful dephasing for {:.2f} steps",
      // step);
    }

    for (long i = 0; i < dephaseBufferLength; i++) {
      delete dephaseBuffer[i];
    }

    if ((params->parrepDephaseLoopStop) &&
        (loop > params->parrepDephaseLoopMax)) {
      SPDLOG_LOGGER_DEBUG(
          log,
          "Reach dephase loop maximum, stop dephasing! Dephased for {} steps",
          step);
      break;
    }
    SPDLOG_LOGGER_DEBUG(log, "Successfully Dephased for {:.2f} fs",
                        step * params->mdTimeStep * params->timeUnit);
  }
}

bool SafeHyperJob::checkState(Matter *current, Matter *reactant) {
  Matter tmp(*current);
  tmp.relax(true);
  if (tmp.compare(*reactant)) {
    return false;
  }
  return true;
}

long SafeHyperJob::refine(Matter *buff[], long length, Matter *reactant) {
  // SPDLOG_LOGGER_DEBUG(log, "[Parallel Replica] Refining transition
  // time.\n");

  bool midTest;
  long min, max, mid;

  min = 0;
  max = length - 1;

  while ((max - min) > 1) {

    mid = min + (max - min) / 2;
    midTest = checkState(buff[mid], reactant);

    if (midTest == false) {
      min = mid;
    } else if (midTest == true) {
      max = mid;
    } else {
      SPDLOG_LOGGER_CRITICAL(log, "Refine step failed!");
      std::exit(1);
    }
  }

  return (min + max) / 2 + 1;
}
