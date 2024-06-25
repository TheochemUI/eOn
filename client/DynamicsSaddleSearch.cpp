#include "DynamicsSaddleSearch.h"
#include "BondBoost.h"
#include "Dimer.h"
#include "Dynamics.h"
#include "ImprovedDimer.h"
#include "Lanczos.h"
#include "LowestEigenmode.h"
#include "MinModeSaddleSearch.h"
#include "NudgedElasticBand.h"

int DynamicsSaddleSearch::run(void) {
  std::vector<std::shared_ptr<Matter>> MDSnapshots;
  std::vector<double> MDTimes;
  SPDLOG_LOGGER_DEBUG(log, "Starting dynamics NEB saddle search");

  ifstream massFile("masses.dat");
  if (massFile.is_open()) {
    SPDLOG_LOGGER_DEBUG(log, "Found mass weights file");
    massFile.close();
    Eigen::VectorXd masses =
        helper_functions::loadMasses("masses.dat", saddle->numberOfAtoms());
    saddle->setMasses(masses);
    SPDLOG_LOGGER_DEBUG(log, "Applied mass weights");
  } else {
    SPDLOG_LOGGER_DEBUG(log, "No mass weights file found");
    massFile.close();
  }

  Dynamics dyn(saddle.get(), params.get());
  SPDLOG_LOGGER_DEBUG(
      log, "Initializing velocities from Maxwell-Boltzmann distribution");
  dyn.setTemperature(params->saddleDynamicsTemperature);
  dyn.setThermalVelocity();

  int dephaseSteps =
      int(floor(params->parrepDephaseTime / params->mdTimeStep + 0.5));

  while (true) {

    SPDLOG_LOGGER_DEBUG(log, "Dephasing: {} steps", dephaseSteps);
    // always start from the initial configuration
    *saddle = *reactant;
    dyn.setThermalVelocity();

    // Dephase MD trajectory
    for (int step = 1; step <= dephaseSteps; step++) {
      dyn.oneStep(step);
    }

    // Check to see if a transition occured
    Matter min(pot, params);
    min = *saddle;
    min.relax();

    if (min.compare(*reactant)) {
      SPDLOG_LOGGER_DEBUG(log, "Dephasing successful");
      break;
    } else {
      SPDLOG_LOGGER_DEBUG(log,
                          "Transition occured during dephasing; Restarting");
      dephaseSteps /= 2;
      if (dephaseSteps < 1)
        dephaseSteps = 1;
    }
  }

  BondBoost bondBoost(saddle.get(), params.get());
  if (params->biasPotential == Hyperdynamics::BOND_BOOST) {
    SPDLOG_LOGGER_DEBUG(log, "Initializing Bond Boost");
    bondBoost.initialize();
  }

  int checkInterval =
      int(params->saddleDynamicsStateCheckInterval / params->mdTimeStep + 0.5);
  int recordInterval =
      int(params->saddleDynamicsRecordInterval / params->mdTimeStep + 0.5);

  if (params->writeMovies == true) {
    saddle->matter2con("dynamics", false);
  }

  for (int step = 1; step <= params->mdSteps; step++) {
    dyn.oneStep(step);

    if (step % recordInterval == 0 && recordInterval != 0) {
      SPDLOG_LOGGER_DEBUG(log, "recording configuration at step {} time {:.3f}",
                          step, step * params->mdTimeStep * params->timeUnit);
      auto tmp = std::shared_ptr<Matter>(saddle);
      *tmp = *saddle;
      MDSnapshots.push_back(tmp);
      MDTimes.push_back(step * params->mdTimeStep);
    }

    if (params->writeMovies == true) {
      saddle->matter2con("dynamics", true);
    }

    if (step % checkInterval == 0) {
      SPDLOG_LOGGER_DEBUG(log, "Minimizing trajectory, step {}", step);

      product = saddle;
      product->relax(false, false);

      if (!product->compare(*reactant)) {
        // SPDLOG_LOGGER_DEBUG(log, "Force calls total: %i",
        // Potential::fcallsTotal);
        SPDLOG_LOGGER_DEBUG(log, "Found new state");
        int image = refineTransition(MDSnapshots, product);
        *saddle = *MDSnapshots[image];
        SPDLOG_LOGGER_DEBUG(log, "Found transition at snapshot image {}",
                            image);
        for (int ii = 0; ii < (int)MDTimes.size(); ii++)
          SPDLOG_LOGGER_DEBUG(log, "MDTimes[{}] = {:.3f}", ii,
                              MDTimes[ii] * params->timeUnit);
        // subtract off half the record interval in order to not introduce a
        // systematic bias towards longer times.
        time = MDTimes[image] - params->saddleDynamicsRecordInterval / 2.0;
        SPDLOG_LOGGER_DEBUG(log, "Transition time {:.2f} fs",
                            time * params->timeUnit);

        NudgedElasticBand neb(reactant, product, params, pot);

        if (params->saddleDynamicsLinearInterpolation == false) {
          SPDLOG_LOGGER_DEBUG(
              log, "Interpolating initial band through MD transition state");
          AtomMatrix reactantToSaddle =
              saddle->pbc(saddle->getPositions() - reactant->getPositions());
          AtomMatrix saddleToProduct =
              saddle->pbc(product->getPositions() - saddle->getPositions());
          SPDLOG_LOGGER_DEBUG(log,
                              "Initial band saved to neb_initial_band.con");
          neb.path[0]->matter2con("neb_initial_band.con", false);
          for (int image = 1; image <= neb.numImages; image++) {
            int mid = neb.numImages / 2 + 1;
            if (image < mid) {
              double frac = ((double)image) / ((double)mid);
              neb.path[image]->setPositions(reactant->getPositions() +
                                            frac * reactantToSaddle);
            } else if (image > mid) {
              double frac =
                  (double)(image - mid) / (double)(neb.numImages - mid + 1);
              neb.path[image]->setPositions(saddle->getPositions() +
                                            frac * saddleToProduct);
            } else if (image == mid) {
              neb.path[image]->setPositions(saddle->getPositions());
            }
            neb.path[image]->matter2con("neb_initial_band.con", true);
          }
          neb.path[neb.numImages + 1]->matter2con("neb_initial_band.con", true);
        } else {
          SPDLOG_LOGGER_DEBUG(
              log, "Linear interpolation between minima used for initial band");
          neb.path[0]->matter2con("neb_initial_band.con", false);
          for (int j = 1; j <= neb.numImages + 1; j++) {
            neb.path[j]->matter2con("neb_initial_band.con", true);
          }
        }

        AtomMatrix mode;
        if (params->nebMaxIterations > 0) {
          LowestEigenmode *minModeMethod;
          if (params->saddleMinmodeMethod == LowestEigenmode::MINMODE_DIMER) {
            if (params->dimerImproved) {
              minModeMethod = new ImprovedDimer(saddle, params, pot);
            } else {
              minModeMethod = new Dimer(saddle, params, pot);
            }
          } else if (params->saddleMinmodeMethod ==
                     LowestEigenmode::MINMODE_LANCZOS) {
            minModeMethod = new Lanczos(saddle, params, pot);
          }

          neb.compute();
          neb.printImageData(true);
          int extremumImage = -1;
          int j;
          for (j = 0; j < neb.numExtrema; j++) {
            //                        if (neb.extremumCurvature[j] < 0.0) {
            if (neb.extremumCurvature[j] <
                params->saddleDynamicsMaxInitCurvature) {
              extremumImage = (int)floor(neb.extremumPosition[j]);
              *saddle = *neb.path[extremumImage];
              double interpDistance =
                  neb.extremumPosition[j] - (double)extremumImage;
              AtomMatrix bandDirection =
                  saddle->pbc(neb.path[extremumImage + 1]->getPositions() -
                              neb.path[extremumImage]->getPositions());
              saddle->setPositions(interpDistance * bandDirection +
                                   saddle->getPositions());
              mode = saddle->pbc(neb.path[extremumImage + 1]->getPositions() -
                                 saddle->getPositions());
              mode.normalize();
              minModeMethod->compute(saddle, mode);
              double eigenvalue = minModeMethod->getEigenvalue();
              SPDLOG_LOGGER_DEBUG(log, "extrema #%i has eigenvalue %.8f", j + 1,
                                  eigenvalue);

              if (eigenvalue < 0) {
                SPDLOG_LOGGER_DEBUG(
                    log, "chose image %i (extrema #%i) as extremum image",
                    extremumImage, j + 1);
                break;
              } else {
                extremumImage = -1;
              }
            }
          }

          delete minModeMethod;

          if (extremumImage != -1) {
            *saddle = *neb.path[extremumImage];
            double interpDistance =
                neb.extremumPosition[j] - (double)extremumImage;
            SPDLOG_LOGGER_DEBUG(log, "interpDistance %f", interpDistance);
            AtomMatrix bandDirection =
                saddle->pbc(neb.path[extremumImage + 1]->getPositions() -
                            neb.path[extremumImage]->getPositions());
            saddle->setPositions(interpDistance * bandDirection +
                                 saddle->getPositions());
            mode = saddle->pbc(neb.path[extremumImage + 1]->getPositions() -
                               saddle->getPositions());
            mode.normalize();
          } else {
            SPDLOG_LOGGER_DEBUG(
                log, "no maxima found, using max energy non-endpoint image");
            double maxEnergy = -INFINITY;
            for (int image = 1; image <= neb.numImages; image++) {
              double U = neb.path[image]->getPotentialEnergy();
              if (U > maxEnergy) {
                maxEnergy = U;
                *saddle = *neb.path[image];
                mode = saddle->pbc(neb.path[image + 1]->getPositions() -
                                   saddle->getPositions());
                mode.normalize();
              }
            }
            if (maxEnergy <= reactant->getPotentialEnergy()) {
              SPDLOG_LOGGER_DEBUG(log, "warning: no barrier found");
              return MinModeSaddleSearch::STATUS_BAD_NO_BARRIER;
            }
          }
        } else {
          neb.maxEnergyImage = neb.numImages / 2 + 1;
        }

        SPDLOG_LOGGER_DEBUG(
            log, "Initial saddle guess saved to saddle_initial_guess.con");
        saddle->matter2con("saddle_initial_guess.con");
        MinModeSaddleSearch search = MinModeSaddleSearch(
            saddle, mode, reactant->getPotentialEnergy(), params, pot);
        int minModeStatus = search.run();

        if (minModeStatus != MinModeSaddleSearch::STATUS_GOOD) {
          SPDLOG_LOGGER_DEBUG(log, "error in min mode saddle search");
          return minModeStatus;
        }

        eigenvalue = search.getEigenvalue();
        eigenvector = search.getEigenvector();
        SPDLOG_LOGGER_DEBUG(log, "eigenvalue: {:.3f}", eigenvalue);

        double barrier =
            saddle->getPotentialEnergy() - reactant->getPotentialEnergy();
        SPDLOG_LOGGER_DEBUG(log, "found barrier of {:.3f}", barrier);
        MDSnapshots.clear();
        MDTimes.clear();
        // SPDLOG_LOGGER_DEBUG(log, "Force calls total: %i",
        // Potential::fcallsTotal);
        return MinModeSaddleSearch::STATUS_GOOD;
      } else {
        SPDLOG_LOGGER_DEBUG(log, "Still in original state");
        MDTimes.clear();
        MDSnapshots.clear();
      }
    }
  }

  MDSnapshots.clear();
  time = params->mdSteps * params->mdTimeStep;
  return MinModeSaddleSearch::STATUS_BAD_MD_TRAJECTORY_TOO_SHORT;
}

int DynamicsSaddleSearch::refineTransition(
    std::vector<std::shared_ptr<Matter>> MDSnapshots,
    std::shared_ptr<Matter> product) {
  int min, max, mid;
  bool midTest;
  min = 0;
  max = MDSnapshots.size() - 1;
  if (max == 0) {
    return 0;
  }

  SPDLOG_LOGGER_DEBUG(log, "refining transition time");

  while ((max - min) > 1) {
    mid = min + (max - min) / 2;
    SPDLOG_LOGGER_DEBUG(log, "minimizing image {}", mid);
    Matter snapshot(pot, params);
    snapshot = *MDSnapshots[mid];

    snapshot.relax(false);

    midTest = snapshot.compare(*reactant);

    if (midTest) {
      SPDLOG_LOGGER_DEBUG(log, "image {} minimizes to reactant", mid);
      min = mid;
    } else {
      SPDLOG_LOGGER_DEBUG(log, "image {} minimizes to product", mid);
      *product = snapshot;
      max = mid;
    }
  }

  return (min + max) / 2;
}

double DynamicsSaddleSearch::getEigenvalue() { return eigenvalue; }

AtomMatrix DynamicsSaddleSearch::getEigenvector() { return eigenvector; }
