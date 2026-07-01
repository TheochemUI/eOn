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

#include "EonLogger.h"
#include <cmath>
#include <format>
#include <fstream>
#include <stdexcept>
#include <string>

#include "BasinHoppingJob.h"
#include "Dynamics.h"
#include "HelperFunctions.h"
#include "ObjectiveFunction.h"
#include "Optimizer.h"
#include "Potential.h"

using namespace eonc::helpers;

std::vector<std::string> BasinHoppingJob::run() {
  bool swapMove;
  double swap_accept = 0.0;
  jump_count = 0; // count of jump movies
  swap_count = 0; // count of swap moves
  disp_count = 0; // count of displacement moves
  int consecutive_rejected_trials = 0;
  double totalAccept = 0.0;
  std::unique_ptr<Matter> minTrial = std::make_unique<Matter>(pot, params);
  std::unique_ptr<Matter> swapTrial = std::make_unique<Matter>(pot, params);

  std::string conFilename = getRelevantFile(params.main_options.conFilename);
  if (!eonc::io::io_ok(current->con2matter(conFilename))) {
    QUILL_LOG_CRITICAL(log, "Failed to load {}", conFilename);
    throw std::runtime_error("failed to load " + conFilename);
  }

  // Sanity Check
  std::vector<long> Elements;
  Elements = getElements(current.get());
  if (params.basin_hopping_options.swap_probability > 0 &&
      Elements.size() == 1) {
    log = eonc::log::traceback();
    QUILL_LOG_CRITICAL(log,
                       "error: [Basin Hopping] swap move probability must be "
                       "zero if there is only one element type\n");
    std::exit(1);
  }

  double randomProb =
      params.basin_hopping_options.initial_random_structure_probability;
  if (randomProb > 0.0) {
    QUILL_LOG_DEBUG(log, "generating random structure with probability {:.4f}",
                    randomProb);
  }
  double u = eonc::helpers::random();
  if (u < params.basin_hopping_options.initial_random_structure_probability) {
    AtomMatrix randomPositions = current->getPositionsFree();
    for (int i = 0; i < current->numberOfFreeAtoms(); i++) {
      for (int j = 0; j < 3; j++) {
        randomPositions(i, j) = eonc::helpers::random();
      }
    }
    randomPositions *= current->getCell();
    current->setPositionsFree(randomPositions);

    pushApart(current, params.basin_hopping_options.push_apart_distance);
  }

  *trial = *current;
  *minTrial = *current;

  current->relax(true);

  double currentEnergy = current->getPotentialEnergy();
  double minimumEnergy = currentEnergy;

  auto minimumEnergyStructure = std::make_shared<Matter>(pot, params);
  *minimumEnergyStructure = *current;
  int nsteps = params.basin_hopping_options.steps +
               params.basin_hopping_options.quenching_steps;
  long totalfc;
  std::ofstream bhFile("bh.dat");

  QUILL_LOG_DEBUG(
      log, "[Basin Hopping] {:4s} {:12s} {:12s} {:12s} {:4s} {:5s} {:5s}",
      "step", "current", "trial", "global min", "fc", "ar", "md");
  QUILL_LOG_DEBUG(
      log, "[Basin Hopping] {:4s} {:12s} {:12s} {:12s} {:4s} {:5s} {:5s}",
      "----", "-------", "-----", "----------", "--", "--", "--");

  int recentAccept = 0;
  double curDisplacement = params.basin_hopping_options.displacement;

  for (int step = 0; step < nsteps; step++) {

    // Swap or displace
    if (randomDouble(1.0) < params.basin_hopping_options.swap_probability &&
        step < params.basin_hopping_options.steps) {
      *swapTrial = *current;
      randomSwap(swapTrial.get());
      swapMove = true;
      *minTrial = *swapTrial;
    } else {
      AtomMatrix displacement;
      displacement = displaceRandom(curDisplacement);

      trial->setPositions(current->getPositions() + displacement);
      swapMove = false;
      pushApart(trial, params.basin_hopping_options.push_apart_distance);

      *minTrial = *trial;
    }

    if (params.debug_options.write_movies) {
      if (!eonc::io::io_ok(trial->matter2con("trials", true))) {
        QUILL_LOG_WARNING(log, "Failed to append trials movie frame");
      }
    }

    // Potential::fcalls = 0;
    minTrial->relax(true);
    // int minfcalls = Potential::fcalls;

    double deltaE = minTrial->getPotentialEnergy() - currentEnergy;
    double p = 0.0;
    if (step >= params.basin_hopping_options.steps) {
      if (deltaE <= 0.0) {
        p = 1.0;
      }
    } else {
      if (deltaE <= 0.0) {
        p = 1.0;
      } else {
        p = std::exp(-deltaE /
                     (params.main_options.temperature * 8.6173324e-5));
      }
    }

    bool accepted = false;
    if (randomDouble(1.0) < p) {
      accepted = true;
      if (params.basin_hopping_options.significant_structure) {
        *current = *minTrial;
      } else {
        *current = *trial;
      }
      if (swapMove) {
        swap_accept += 1;
      }
      if (step < params.basin_hopping_options.steps) {
        totalAccept += 1;
        recentAccept += 1;
      }

      currentEnergy = minTrial->getPotentialEnergy();

      if (currentEnergy < minimumEnergy) {
        minimumEnergy = currentEnergy;
        *minimumEnergyStructure = *minTrial;
        if (!eonc::io::io_ok(minimumEnergyStructure->matter2con("min.con"))) {
          QUILL_LOG_WARNING(log, "Failed to write min.con");
        }
      }

      if (params.basin_hopping_options.write_unique) {
        bool newStructure = true;
        for (unsigned int i = 0; i < uniqueEnergies.size(); i++) {
          // if minTrial has a different energy or a different structure
          // it is new, otherwise it is old
          if (std::fabs(currentEnergy - uniqueEnergies[i]) <
              params.structure_comparison_options.energy_difference) {
            if (current->compare(*uniqueStructures[i],
                                 params.structure_comparison_options
                                     .indistinguishable_atoms)) {
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
          if (!eonc::io::io_ok(current->matter2con(fname))) {
            QUILL_LOG_WARNING(log, "Failed to write {}", fname);
          }
          returnFiles.push_back(fname);

          snprintf(fname, 128, "energy_%.5i.dat", step + 1);
          returnFiles.push_back(fname);
          {
            std::ofstream fh(fname);
            if (fh)
              fh << std::format("{:.10e}\n", currentEnergy);
          }
        }
      }

      consecutive_rejected_trials = 0; // STC: I think this should go here.
    } else {
      consecutive_rejected_trials++;
    }

    if (params.debug_options.write_movies) {
      if (!eonc::io::io_ok(minTrial->matter2con("movie", true))) {
        QUILL_LOG_WARNING(log, "Failed to append basin-hopping movie frame");
      }
    }

    // totalfc = Potential::fcallsTotal;
    char acceptReject[2];
    acceptReject[1] = '\0';
    if (accepted) {
      acceptReject[0] = 'A';
    } else {
      acceptReject[0] = 'R';
    }
    // QUILL_LOG_DEBUG(log, "[Basin Hopping] %5i %12.3f %12.3f %12.3f %4i
    // %5.3f %5.3f %1s\n",
    //        step+1, currentEnergy, minTrial->getPotentialEnergy(),
    //        minimumEnergy, minfcalls, totalAccept/((double)step+1),
    //        curDisplacement, acceptReject);
    // fprintf(pFile, "%6i %9ld %12.4e %12.4e\n",step+1,totalfc,currentEnergy,
    // minTrial->getPotentialEnergy());

    if (minimumEnergy < params.basin_hopping_options.stop_energy) {
      break;
    }

    if (consecutive_rejected_trials == params.basin_hopping_options.jump_max &&
        step < params.basin_hopping_options.steps) {
      consecutive_rejected_trials = 0;
      AtomMatrix jump;
      for (int j = 0; j < params.basin_hopping_options.jump_steps; j++) {
        jump_count++;
        jump = displaceRandom(curDisplacement);
        current->setPositions(current->getPositions() + jump);
        if (params.basin_hopping_options.significant_structure) {
          pushApart(current, params.basin_hopping_options.push_apart_distance);
          current->relax(true);
        }
        currentEnergy = current->getPotentialEnergy();
        if (currentEnergy < minimumEnergy) {
          minimumEnergy = currentEnergy;
          *minimumEnergyStructure = *current;
        }
      }
    }

    int nadjust = params.basin_hopping_options.adjust_period;
    double adjustFraction = params.basin_hopping_options.adjust_fraction;
    if ((step + 1) % nadjust == 0 &&
        params.basin_hopping_options.adjust_displacement) {
      double recentRatio =
          static_cast<double>(recentAccept) / static_cast<double>(nadjust);
      if (recentRatio > params.basin_hopping_options.target_ratio) {
        curDisplacement *= 1.0 + adjustFraction;
      } else {
        curDisplacement *= 1.0 - adjustFraction;
      }

      // QUILL_LOG_DEBUG(log, "recentRatio %.3f md: %.3f\n", recentRatio,
      // curDisplacement);
      recentAccept = 0;
    }
  }
  bhFile.close();

  /* Save Results */

  std::string resultsFilename("results.dat");
  returnFiles.push_back(resultsFilename);

  if (params.debug_options.write_movies) {
    std::string movieFilename("movie.xyz");
    returnFiles.push_back(movieFilename);
  }

  {
    std::ofstream out(resultsFilename, std::ios::binary);
    if (out) {
      out << std::format("{} termination_reason\n", 0);
      out << "GOOD termination_reason_text\n";
      out << std::format("{:.6f} minimum_energy\n", minimumEnergy);
      out << std::format("{} random_seed\n", params.main_options.randomSeed);
      out << std::format("{:.3f} acceptance_ratio\n",
                         totalAccept / params.basin_hopping_options.steps);
      if (params.basin_hopping_options.swap_probability > 0) {
        out << std::format("{:.3f} swap_acceptance_ratio\n",
                           swap_accept / static_cast<double>(swap_count));
      }
      out << std::format("{} total_normal_displacement_steps\n",
                         disp_count - jump_count -
                             params.basin_hopping_options.quenching_steps);
      out << std::format("{} total_jump_steps\n", jump_count);
      out << std::format("{} total_swap_steps\n", swap_count);
      out << std::format("{} total_force_calls\n",
                         PotRegistry::get().total_force_calls());
    }
  }

  std::string productFilename("min.con");
  returnFiles.push_back(productFilename);
  if (!eonc::io::io_ok(minimumEnergyStructure->matter2con(productFilename))) {
    QUILL_LOG_ERROR(log, "Failed to write {}", productFilename);
  }

  std::string bhFilename("bh.dat");
  returnFiles.push_back(bhFilename);

  // minTrial and swapTrial automatically cleaned up by unique_ptr
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
  if (params.basin_hopping_options.single_atom_displace) {
    m = randomInt(0, trial->numberOfAtoms() - 1);
    num = m + 1;
  }

  for (int i = m; i < num; i++) {
    double dist = distvec(i);
    double disp = 0.0; // displacement size, possibly scaled

    if (!trial->getFixed(i)) {
      if (params.basin_hopping_options.displacement_algorithm == "standard") {
        disp = curDisplacement;
      }
      // scale displacement linearly with the particle radius
      else if (params.basin_hopping_options.displacement_algorithm ==
               "linear") {
        double Cs = curDisplacement / distvec.maxCoeff();
        disp = Cs * dist;
      }
      // scale displacement quadratically with the particle radius
      else if (params.basin_hopping_options.displacement_algorithm ==
               "quadratic") {
        double Cq = curDisplacement / (distvec.maxCoeff() * distvec.maxCoeff());
        disp = Cq * dist * dist;
      } else {
        log = eonc::log::traceback();
        QUILL_LOG_CRITICAL(log, "Unknown displacement_algorithm\n");
        std::exit(1);
      }
      for (int j = 0; j < 3; j++) {
        if (params.basin_hopping_options.displacement_distribution ==
            "uniform") {
          displacement(i, j) = randomDouble(2 * disp) - disp;
        } else if (params.basin_hopping_options.displacement_distribution ==
                   "gaussian") {
          displacement(i, j) = gaussRandom(0.0, disp);
        } else {
          log = eonc::log::traceback();
          QUILL_LOG_CRITICAL(log, "Unknown displacement_distribution\n");
          std::exit(1);
        }
      }
    }
  }
  return displacement;
}

void BasinHoppingJob::randomSwap(Matter *matter) {
  swap_count++;
  std::vector<long> Elements;
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

std::vector<long> BasinHoppingJob::getElements(Matter *matter) {
  std::array<int, 118> allElements{};
  std::vector<long> elements;

  for (long y = 0; y < matter->numberOfAtoms(); ++y) {
    if (!matter->getFixed(y)) {
      const int index = matter->getAtomicNr(y);
      allElements[index] = 1;
    }
  }

  for (int i = 0; i < 118; ++i) {
    if (allElements[i] != 0) {
      elements.push_back(i);
    }
  }

  return elements;
}

VectorXd BasinHoppingJob::calculateDistanceFromCenter(Matter *matter) {
  AtomMatrix pos = matter->getPositions();
  Vector3d cen(0, 0, 0);
  int num = matter->numberOfAtoms();

  cen = pos.colwise().sum() / static_cast<double>(num);

  VectorXd dist(num);

  for (int n = 0; n < num; n++) {
    pos.row(n) -= cen;
    dist(n) = pos.row(n).norm();
  }

  return dist;
}
