/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
**
** Copyright (c) 2010--present, eOn Development Team
** All rights reserved.
**
** Repo:
** https://github.com/TheochemUI/eOn
*/

#include <math.h>
#include <stdio.h>
#include <string>

#include "BasinHoppingJob.h"
#include "Dynamics.h"
#include "HelperFunctions.h"
#include "ObjectiveFunction.h"
#include "Optimizer.h"
#include "Potential.h"

using namespace std;
using namespace helper_functions;

std::vector<std::string> BasinHoppingJob::run(void) {
  bool swapMove;
  double swap_accept = 0.0;
  jump_count = 0; // count of jump movies
  swap_count = 0; // count of swap moves
  disp_count = 0; // count of displacement moves
  int consecutive_rejected_trials = 0;
  double totalAccept = 0.0;
  Matter *minTrial = new Matter(pot, params);
  Matter *swapTrial = new Matter(pot, params);

  string conFilename = getRelevantFile(params->main_options.conFilename);
  current->con2matter(conFilename);

  // Sanity Check
  vector<long> Elements;
  Elements = getElements(current.get());
  if (params->basin_hopping_options.swap_probability > 0 && Elements.size() == 1) {
    log = spdlog::get("_traceback");
    SPDLOG_LOGGER_CRITICAL(
        log, "error: [Basin Hopping] swap move probability must be "
             "zero if there is only one element type\n");
    std::exit(1);
  }

  double randomProb = params->basin_hopping_options.initial_random_structure_probability;
  if (randomProb > 0.0) {
    SPDLOG_LOGGER_DEBUG(
        log, "generating random structure with probability {:.4f}", randomProb);
  }
  double u = helper_functions::random();
  if (u < params->basin_hopping_options.initial_random_structure_probability) {
    AtomMatrix randomPositions = current->getPositionsFree();
    for (int i = 0; i < current->numberOfFreeAtoms(); i++) {
      for (int j = 0; j < 3; j++) {
        randomPositions(i, j) = helper_functions::random();
      }
    }
    randomPositions *= current->getCell();
    current->setPositionsFree(randomPositions);

    pushApart(current, params->basin_hopping_options.push_apart_distance);
  }

  *trial = *current;
  *minTrial = *current;

  current->relax(true);

  double currentEnergy = current->getPotentialEnergy();
  double minimumEnergy = currentEnergy;

  auto minimumEnergyStructure = std::make_shared<Matter>(pot, params);
  *minimumEnergyStructure = *current;
  int nsteps = params->basin_hopping_options.steps + params->basin_hopping_options.quenching_steps;
  long totalfc;
  FILE *pFile;
  pFile = fopen("bh.dat", "w");

  SPDLOG_LOGGER_DEBUG(
      log, "[Basin Hopping] {:4s} {:12s} {:12s} {:12s} {:4s} {:5s} {:5s}",
      "step", "current", "trial", "global min", "fc", "ar", "md");
  SPDLOG_LOGGER_DEBUG(
      log, "[Basin Hopping] {:4s} {:12s} {:12s} {:12s} {:4s} {:5s} {:5s}",
      "----", "-------", "-----", "----------", "--", "--", "--");

  int recentAccept = 0;
  double curDisplacement = params->basin_hopping_options.displacement;

  for (int step = 0; step < nsteps; step++) {

    // Swap or displace
    if (randomDouble(1.0) < params->basin_hopping_options.swap_probability &&
        step < params->basin_hopping_options.steps) {
      *swapTrial = *current;
      randomSwap(swapTrial);
      swapMove = true;
      *minTrial = *swapTrial;
    } else {
      AtomMatrix displacement;
      displacement = displaceRandom(curDisplacement);

      trial->setPositions(current->getPositions() + displacement);
      swapMove = false;
      pushApart(trial, params->basin_hopping_options.push_apart_distance);

      *minTrial = *trial;
    }

    if (params->debug_options.write_movies == true) {
      trial->matter2con("trials", true);
    }

    // Potential::fcalls = 0;
    minTrial->relax(true);
    // int minfcalls = Potential::fcalls;

    double deltaE = minTrial->getPotentialEnergy() - currentEnergy;
    double p = 0.0;
    if (step >= params->basin_hopping_options.steps) {
      if (deltaE <= 0.0) {
        p = 1.0;
      }
    } else {
      if (deltaE <= 0.0) {
        p = 1.0;
      } else {
        p = exp(-deltaE / (params->main_options.temperature * 8.6173324e-5));
      }
    }

    bool accepted = false;
    if (randomDouble(1.0) < p) {
      accepted = true;
      if (params->basin_hopping_options.significant_structure) {
        *current = *minTrial;
      } else {
        *current = *trial;
      }
      if (swapMove) {
        swap_accept += 1;
      }
      if (step < params->basin_hopping_options.steps) {
        totalAccept += 1;
        recentAccept += 1;
      }

      currentEnergy = minTrial->getPotentialEnergy();

      if (currentEnergy < minimumEnergy) {
        minimumEnergy = currentEnergy;
        *minimumEnergyStructure = *minTrial;
        minimumEnergyStructure->matter2con("min.con");
      }

      if (params->basin_hopping_options.write_unique) {
        bool newStructure = true;
        for (unsigned int i = 0; i < uniqueEnergies.size(); i++) {
          // if minTrial has a different energy or a different structure
          // it is new, otherwise it is old
          if (fabs(currentEnergy - uniqueEnergies[i]) <
              params->structure_comparison_options.energy_difference) {
            if (current->compare(*uniqueStructures[i],
                                 params->structure_comparison_options.indistinguishable_atoms) == true) {
              newStructure = false;
            }
          }
        }

        if (newStructure) {
          uniqueEnergies.push_back(currentEnergy);
          auto currentCopy = std::make_shared<Matter>(pot, params);
          *currentCopy = *current;
          uniqueStructures.push_back(currentCopy);

          char fname[128];
          snprintf(fname, 128, "min_%.5i.con", step + 1);
          current->matter2con(fname);
          returnFiles.push_back(fname);

          snprintf(fname, 128, "energy_%.5i.dat", step + 1);
          returnFiles.push_back(fname);
          FILE *fh = fopen(fname, "w");
          fprintf(fh, "%.10e\n", currentEnergy);
          fclose(fh);
        }
      }

      consecutive_rejected_trials = 0; // STC: I think this should go here.
    } else {
      consecutive_rejected_trials++;
    }

    if (params->debug_options.write_movies == true) {
      minTrial->matter2con("movie", true);
    }

    // totalfc = Potential::fcallsTotal;
    char acceptReject[2];
    acceptReject[1] = '\0';
    if (accepted) {
      acceptReject[0] = 'A';
    } else {
      acceptReject[0] = 'R';
    }
    // SPDLOG_LOGGER_DEBUG(log, "[Basin Hopping] %5i %12.3f %12.3f %12.3f %4i
    // %5.3f %5.3f %1s\n",
    //        step+1, currentEnergy, minTrial->getPotentialEnergy(),
    //        minimumEnergy, minfcalls, totalAccept/((double)step+1),
    //        curDisplacement, acceptReject);
    // fprintf(pFile, "%6i %9ld %12.4e %12.4e\n",step+1,totalfc,currentEnergy,
    // minTrial->getPotentialEnergy());

    if (minimumEnergy < params->basin_hopping_options.stop_energy) {
      break;
    }

    if (consecutive_rejected_trials == params->basin_hopping_options.jump_max &&
        step < params->basin_hopping_options.steps) {
      consecutive_rejected_trials = 0;
      AtomMatrix jump;
      for (int j = 0; j < params->basin_hopping_options.jump_steps; j++) {
        jump_count++;
        jump = displaceRandom(curDisplacement);
        current->setPositions(current->getPositions() + jump);
        if (params->basin_hopping_options.significant_structure) {
          pushApart(current, params->basin_hopping_options.push_apart_distance);
          current->relax(true);
        }
        currentEnergy = current->getPotentialEnergy();
        if (currentEnergy < minimumEnergy) {
          minimumEnergy = currentEnergy;
          *minimumEnergyStructure = *current;
        }
      }
    }

    int nadjust = params->basin_hopping_options.adjust_period;
    double adjustFraction = params->basin_hopping_options.adjust_fraction;
    if ((step + 1) % nadjust == 0 &&
        params->basin_hopping_options.adjust_displacement == true) {
      double recentRatio = ((double)recentAccept) / ((double)nadjust);
      if (recentRatio > params->basin_hopping_options.target_ratio) {
        curDisplacement *= 1.0 + adjustFraction;
      } else {
        curDisplacement *= 1.0 - adjustFraction;
      }

      // SPDLOG_LOGGER_DEBUG(log, "recentRatio %.3f md: %.3f\n", recentRatio,
      // curDisplacement);
      recentAccept = 0;
    }
  }
  fclose(pFile);

  /* Save Results */

  FILE *fileResults, *fileProduct;

  std::string resultsFilename("results.dat");
  returnFiles.push_back(resultsFilename);
  fileResults = fopen(resultsFilename.c_str(), "wb");

  if (params->debug_options.write_movies == true) {
    std::string movieFilename("movie.xyz");
    returnFiles.push_back(movieFilename);
  }

  fprintf(fileResults, "%d termination_reason\n", 0);
  fprintf(fileResults, "%.6f minimum_energy\n", minimumEnergy);
  fprintf(fileResults, "%ld random_seed\n", params->main_options.randomSeed);
  fprintf(fileResults, "%.3f acceptance_ratio\n",
          totalAccept / params->basin_hopping_options.steps);
  if (params->basin_hopping_options.swap_probability > 0) {
    fprintf(fileResults, "%.3f swap_acceptance_ratio\n",
            swap_accept / double(swap_count));
  }
  fprintf(fileResults, "%ld total_normal_displacement_steps\n",
          disp_count - jump_count - params->basin_hopping_options.quenching_steps);
  fprintf(fileResults, "%d total_jump_steps\n", jump_count);
  fprintf(fileResults, "%d total_swap_steps\n", swap_count);
  // fprintf(fileResults, "%d total_force_calls\n", Potential::fcallsTotal);
  fclose(fileResults);

  std::string productFilename("min.con");
  returnFiles.push_back(productFilename);
  fileProduct = fopen(productFilename.c_str(), "wb");
  minimumEnergyStructure->matter2con(fileProduct);
  fclose(fileProduct);

  std::string bhFilename("bh.dat");
  returnFiles.push_back(bhFilename);

  delete minTrial;
  delete swapTrial;
  return returnFiles;
}

AtomMatrix BasinHoppingJob::displaceRandom(double curDisplacement) {
  disp_count++;
  // Create a random displacement
  AtomMatrix displacement;
  displacement.resize(trial->numberOfAtoms(), 3);
  displacement.setZero();
  VectorXd distvec = calculateDistanceFromCenter(current.get());
  int num = trial->numberOfAtoms();
  int m = 0;
  if (params->basin_hopping_options.single_atom_displace) {
    m = randomInt(0, trial->numberOfAtoms() - 1);
    num = m + 1;
  }

  for (int i = m; i < num; i++) {
    double dist = distvec(i);
    double disp = 0.0; // displacement size, possibly scaled

    if (!trial->getFixed(i)) {
      if (params->basin_hopping_options.displacement_algorithm == "standard") {
        disp = curDisplacement;
      }
      // scale displacement linearly with the particle radius
      else if (params->basin_hopping_options.displacement_algorithm == "linear") {
        double Cs = curDisplacement / distvec.maxCoeff();
        disp = Cs * dist;
      }
      // scale displacement quadratically with the particle radius
      else if (params->basin_hopping_options.displacement_algorithm == "quadratic") {
        double Cq = curDisplacement / (distvec.maxCoeff() * distvec.maxCoeff());
        disp = Cq * dist * dist;
      } else {
        log = spdlog::get("_traceback");
        SPDLOG_LOGGER_CRITICAL(log, "Unknown displacement_algorithm\n");
        std::exit(1);
      }
      for (int j = 0; j < 3; j++) {
        if (params->basin_hopping_options.displacement_distribution == "uniform") {
          displacement(i, j) = randomDouble(2 * disp) - disp;
        } else if (params->basin_hopping_options.displacement_distribution == "gaussian") {
          displacement(i, j) = gaussRandom(0.0, disp);
        } else {
          log = spdlog::get("_traceback");
          SPDLOG_LOGGER_CRITICAL(log, "Unknown displacement_distribution\n");
          std::exit(1);
        }
      }
    }
  }
  return displacement;
}

void BasinHoppingJob::randomSwap(Matter *matter) {
  swap_count++;
  vector<long> Elements;
  Elements = getElements(matter);

  long ela;
  long elb;
  long ia = randomInt(0, Elements.size() - 1);
  ela = Elements.at(ia);
  Elements.erase(Elements.begin() + ia);

  long ib = randomInt(0, Elements.size() - 1);
  elb = Elements.at(ib);

  int changera = 0;
  int changerb = 0;

  changera = randomInt(0, matter->numberOfAtoms() - 1);
  while (matter->getAtomicNr(changera) != ela) {
    changera = randomInt(0, matter->numberOfAtoms() - 1);
  }

  changerb = randomInt(0, matter->numberOfAtoms() - 1);
  while (matter->getAtomicNr(changerb) != elb) {
    changerb = randomInt(0, matter->numberOfAtoms() - 1);
  }

  double posax = matter->getPosition(changera, 0);
  double posay = matter->getPosition(changera, 1);
  double posaz = matter->getPosition(changera, 2);

  matter->setPosition(changera, 0, matter->getPosition(changerb, 0));
  matter->setPosition(changera, 1, matter->getPosition(changerb, 1));
  matter->setPosition(changera, 2, matter->getPosition(changerb, 2));

  matter->setPosition(changerb, 0, posax);
  matter->setPosition(changerb, 1, posay);
  matter->setPosition(changerb, 2, posaz);
}

vector<long> BasinHoppingJob::getElements(Matter *matter) {
  int allElements[118] = {0};
  vector<long> Elements;

  for (long y = 0; y < matter->numberOfAtoms(); y++) {
    if (!matter->getFixed(y)) {
      int index = matter->getAtomicNr(y);
      allElements[index] = 1;
    }
  }

  for (int i = 0; i < 118; i++) {
    if (allElements[i] != 0) {
      Elements.push_back(i);
    }
  }

  return Elements;
}

VectorXd BasinHoppingJob::calculateDistanceFromCenter(Matter *matter) {
  AtomMatrix pos = matter->getPositions();
  Vector3d cen(0, 0, 0);
  int num = matter->numberOfAtoms();

  cen = pos.colwise().sum() / (double)num;

  VectorXd dist(num);

  for (int n = 0; n < num; n++) {
    pos.row(n) -= cen;
    dist(n) = pos.row(n).norm();
  }

  return dist;
}
