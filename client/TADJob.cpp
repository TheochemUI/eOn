#include <cstdlib>

#include "BaseStructures.h"
#include "Dynamics.h"
#include "Matter.h"
#include "Optimizer.h"
#include "TADJob.h"

std::vector<std::string> TADJob::run(void) {
  current = std::make_shared<Matter>(pot, params);
  reactant = std::make_shared<Matter>(pot, params);
  saddle = std::make_shared<Matter>(pot, params);
  crossing = std::make_shared<Matter>(pot, params);
  product = std::make_shared<Matter>(pot, params);
  final_state = std::make_shared<Matter>(pot, params);
  final_tmp = std::make_shared<Matter>(pot, params);

  minimizeFCalls = mdFCalls = refineFCalls = dephaseFCalls = 0;
  time = 0.0;
  string reactantFilename =
      helper_functions::getRelevantFile(params->conFilename);
  current->con2matter(reactantFilename);

  SPDLOG_LOGGER_DEBUG(log, "Minimizing initial reactant");
  // long refFCalls = Potential::fcalls;
  *reactant = *current;
  reactant->relax();
  // minimizeFCalls += (Potential::fcalls - refFCalls);

  SPDLOG_LOGGER_DEBUG(log, "Temperature Accelerated Dynamics, running");
  SPDLOG_LOGGER_DEBUG(log,
                      "High temperature MD simulation running at {:.2f} K to "
                      "simulate dynamics at {:.2f} K",
                      params->temperature, params->tadLowT);

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

  return returnFiles;
}

int TADJob::dynamics() {
  bool transitionFlag = false, recordFlag = true, stopFlag = false,
       firstTransitFlag = false;
  long nFreeCoord = reactant->numberOfFreeAtoms() * 3;
  long mdBufferLength, refFCalls;
  long step = 0, refineStep,
       newStateStep = 0; // check that newStateStep is set before used
  long nCheck = 0, nRecord = 0, nState = 0;
  long StateCheckInterval, RecordInterval;
  double kinE, kinT, avgT, varT;
  double kB = params->kB;
  double correctedTime = 0.0;
  double stopTime = 0.0, sumSimulatedTime = 0.0;
  double Temp = 0.0, sumT = 0.0, sumT2 = 0.0;
  double correctionFactor = 1.0;
  double transitionTime_current = 0.0, transitionTime_previous = 0.0;
  double delta, minmu, factor, highT, lowT;

  AtomMatrix velocity;
  AtomMatrix reducedForces;

  minCorrectedTime = 1.0e200;
  lowT = params->tadLowT;
  highT = params->temperature;
  delta = params->tadConfidence;
  minmu = params->tadMinPrefactor;
  factor = std::log(1.0 / delta) / minmu;
  StateCheckInterval =
      int(params->parrepStateCheckInterval / params->mdTimeStep);
  RecordInterval = int(params->parrepRecordInterval / params->mdTimeStep);
  Temp = params->temperature;
  newStateFlag = metaStateFlag = false;

  mdBufferLength = long(StateCheckInterval / RecordInterval);
  std::vector<std::shared_ptr<Matter>> mdBuffer;
  mdBuffer.resize(mdBufferLength);
  for (long i = 0; i < mdBufferLength; i++) {
    mdBuffer[i] = std::make_shared<Matter>(pot, params);
  }
  timeBuffer = new double[mdBufferLength];

  Dynamics TAD(current.get(), params.get());
  TAD.setThermalVelocity();

  // dephase the trajectory so that it is thermal and independent of others
  // refFCalls = Potential::fcalls;
  dephase();
  // dephaseFCalls = Potential::fcalls - refFCalls;

  SPDLOG_LOGGER_DEBUG(
      log,
      "Starting MD run\nTemperature: {:.2f} Kelvin"
      "Total Simulation Time: {:.2f} fs\nTime Step: {:.2f} fs\nTotal Steps: "
      "{}\n",
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

    kinE = current->getKineticEnergy();
    kinT = (2.0 * kinE / nFreeCoord / kB);
    sumT += kinT;
    sumT2 += kinT * kinT;
    // SPDLOG_LOGGER_DEBUG(log, "steps = {:10d} temp = {:10.5f} ", step, kinT);

    TAD.oneStep();
    mdFCalls++;

    time += params->mdTimeStep;
    nCheck++; // count up to params->parrepStateCheckInterval before checking
              // for a transition
    step++;
    // SPDLOG_LOGGER_DEBUG(log, "step = {:4d}, time= {:10.4f}", step, time);
    //  standard conditions; record mater object in the transition buffer
    if (params->parrepRefineTransition && recordFlag && !newStateFlag) {
      if (nCheck % RecordInterval == 0) {
        *mdBuffer[nRecord] = *current;
        timeBuffer[nRecord] = time;
        nRecord++; // current location in the buffer
      }
    }

    // time to do a state check; if a transiton if found, stop recording
    if ((nCheck == StateCheckInterval) && !newStateFlag) {
      nCheck = 0;  // reinitialize check state counter
      nRecord = 0; // restart the buffer
      // refFCalls = Potential::fcalls;
      transitionFlag = checkState(current.get(), reactant.get());
      // minimizeFCalls += Potential::fcalls - refFCalls;
      if (transitionFlag == true) {
        nState++;
        SPDLOG_LOGGER_DEBUG(log, "New State %ld: ", nState);
        *final_tmp = *current;
        transitionTime = time;
        newStateStep = step; // remember the step when we are in a new state
        transitionStep = newStateStep;
        firstTransitFlag = 1;
      }
    }
    // printf("step=%ld, time=%lf, biasPot=%lf",step,time,boostPotential);
    // Refine transition step

    if (transitionFlag) {
      // SPDLOG_LOGGER_DEBUG(log, "[Parallel Replica] Refining transition
      // time.");
      refFCalls = Potential::fcalls;
      refineStep = refine(mdBuffer, mdBufferLength, reactant.get());

      transitionStep =
          newStateStep - StateCheckInterval + refineStep * RecordInterval;
      transitionTime_current = timeBuffer[refineStep];
      *crossing = *mdBuffer[refineStep];
      transitionTime = transitionTime_current - transitionTime_previous;
      transitionTime_previous = transitionTime_current;
      barrier = crossing->getPotentialEnergy() - reactant->getPotentialEnergy();
      SPDLOG_LOGGER_DEBUG(log, "barrier= {:.3f}", barrier);
      correctionFactor = 1.0 * exp(barrier / kB * (1.0 / lowT - 1.0 / highT));
      correctedTime = transitionTime * correctionFactor;
      sumSimulatedTime += transitionTime;

      // reverse the momentum;
      *current = *mdBuffer[refineStep - 1];
      velocity = current->getVelocities();
      velocity = velocity * (-1);
      current->setVelocities(velocity);

      if (correctedTime < minCorrectedTime) {
        minCorrectedTime = correctedTime;
        *saddle = *crossing;
        *final_state = *final_tmp;
      }
      stopTime = factor * pow(minCorrectedTime / factor, lowT / highT);
      SPDLOG_LOGGER_DEBUG(
          log,
          "tranisitonTime= {:.3e} s, Barrier= {:.3f} eV, correctedTime= {:.3e} "
          "s, "
          "SimulatedTime= {:.3e} s, minCorTime= {:.3e} s, stopTime= {:.3e} s",
          transitionTime * 1e-15 * params->timeUnit, barrier,
          correctedTime * 1e-15 * params->timeUnit,
          sumSimulatedTime * 1e-15 * params->timeUnit,
          minCorrectedTime * 1.0e-15 * params->timeUnit,
          stopTime * 1.0e-15 * params->timeUnit);

      // refineFCalls += Potential::fcalls - refFCalls;
      transitionFlag = false;
    }

    // we have run enough md steps; time to stop
    if (firstTransitFlag && sumSimulatedTime >= stopTime) {
      stopFlag = true;
      newStateFlag = true;
    }

    // stdout Progress
    if ((step % tenthSteps == 0) || (step == params->mdSteps)) {
      double maxAtomDistance = current->perAtomNorm(*reactant);
      SPDLOG_LOGGER_DEBUG(
          log, "progress: {:.0f}%, max displacement: {:.3lf}, step {:7d}/{:d}",
          (double)100.0 * step / params->mdSteps, maxAtomDistance, step,
          params->mdSteps);
    }

    if (step == params->mdSteps) {
      stopFlag = true;
      if (firstTransitFlag) {
        SPDLOG_LOGGER_DEBUG(log, "Detected one transition");
      } else {
        SPDLOG_LOGGER_DEBUG(log, "Failed to detect any transition");
      }
    }
  }

  // calculate avearges
  avgT = sumT / step;
  varT = sumT2 / step - avgT * avgT;

  SPDLOG_LOGGER_DEBUG(
      log,
      "Temperature : Average = %lf ; Stddev = %lf ; Factor = %lf; "
      "Average_Boost = %lf",
      avgT, sqrt(varT), varT / avgT / avgT * nFreeCoord / 2,
      minCorrectedTime / step / params->mdTimeStep);
  if (isfinite(avgT) == 0) {
    SPDLOG_LOGGER_DEBUG(log,
                        "Infinite average temperature, something went wrong!");
    newStateFlag = false;
  }

  *product = *final_state;
  // new state was detected; determine refined transition time
  delete[] timeBuffer;

  if (newStateFlag) {
    return 1;
  } else {
    return 0;
  }
}

void TADJob::saveData(int status) {
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
  // fprintf(fileResults, "%ld force_calls_dephase\n", dephaseFCalls);
  fprintf(fileResults, "%ld force_calls_dynamics\n", mdFCalls);
  fprintf(fileResults, "%ld force_calls_minimize\n", minimizeFCalls);
  fprintf(fileResults, "%ld force_calls_refine\n", refineFCalls);

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

void TADJob::dephase() {
  bool transitionFlag = false;
  long step, stepNew, loop;
  long DephaseSteps;
  long dephaseBufferLength, dephaseRefineStep;
  AtomMatrix velocity;

  DephaseSteps = int(params->parrepDephaseTime / params->mdTimeStep);
  Dynamics dephaseDynamics(current.get(), params.get());
  SPDLOG_LOGGER_DEBUG(log, "Dephasing for {:.2f} fs",
                      params->parrepDephaseTime * params->timeUnit);

  step = stepNew = loop = 0;

  while (step < DephaseSteps) {
    // this should be allocated once, and of length DephaseSteps
    dephaseBufferLength = DephaseSteps - step;
    loop++;
    std::vector<std::shared_ptr<Matter>> dephaseBuffer(dephaseBufferLength);

    for (long i = 0; i < dephaseBufferLength; i++) {
      dephaseBuffer[i] = std::make_shared<Matter>(pot, params);
      dephaseDynamics.oneStep();
      *dephaseBuffer[i] = *current;
    }

    transitionFlag = checkState(current.get(), reactant.get());

    if (transitionFlag) {
      dephaseRefineStep =
          refine(dephaseBuffer, dephaseBufferLength, reactant.get());
      SPDLOG_LOGGER_DEBUG(log, "loop = {}; dephase refine step = {}", loop,
                          dephaseRefineStep);
      transitionStep = dephaseRefineStep - 1; // check that this is correct
      transitionStep = (transitionStep > 0) ? transitionStep : 0;
      SPDLOG_LOGGER_DEBUG(
          log,
          "Dephasing warning: in a new state, inverse the momentum and restart "
          "from step %ld",
          step + transitionStep);
      *current = *dephaseBuffer[transitionStep];
      velocity = current->getVelocities();
      velocity = velocity * (-1);
      current->setVelocities(velocity);
      step = step + transitionStep;
    } else {
      step = step + dephaseBufferLength;
      // SPDLOG_LOGGER_DEBUG(log, "Successful dephasing for {.2f} steps ",
      // step);
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

bool TADJob::checkState(Matter *current, Matter *reactant) {
  Matter tmp(pot, params);
  tmp = *current;
  tmp.relax(true);
  if (tmp.compare(*reactant)) {
    return false;
  }
  return true;
}

bool TADJob::saddleSearch(std::shared_ptr<Matter> cross) {
  AtomMatrix mode;
  long status;
  mode = cross->getPositions() - reactant->getPositions();
  mode.normalize();
  dimerSearch = NULL;
  // TODO: Unhandled .get()
  dimerSearch = new MinModeSaddleSearch(
      cross, mode, reactant->getPotentialEnergy(), params, pot);
  status = dimerSearch->run();
  SPDLOG_LOGGER_DEBUG(log, "dimer search status %ld", status);
  if (status != MinModeSaddleSearch::STATUS_GOOD) {
    return false;
  }
  return false;
}

long TADJob::refine(std::vector<std::shared_ptr<Matter>> buff, long length,
                    Matter *reactant) {
  // SPDLOG_LOGGER_DEBUG(log, "[Parallel Replica] Refining transition time.");

  bool midTest;
  long min, max, mid;

  min = 0;
  max = length - 1;

  while ((max - min) > 1) {

    mid = min + (max - min) / 2;
    midTest = checkState(buff[mid].get(), reactant);

    if (midTest == false) {
      min = mid;
    } else if (midTest == true) {
      max = mid;
    } else {
      SPDLOG_LOGGER_DEBUG(log, "Refine step failed!");
      exit(1);
    }
  }

  return (min + max) / 2 + 1;
}
