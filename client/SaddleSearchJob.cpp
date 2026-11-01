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
#include "eon/SaddleSearchJob.h"
#ifdef WITH_ARTN
#include "eon/ARTnSaddleSearch.h"
#endif
#include "eon/EonLogger.h"
#include "eon/EpiCenters.h"
#include "eon/HelperFunctions.h"
#include "eon/Potential.h"
#include "eon/PrfoRefine.h"
#include "eon/XtsciIrc.h"

#include <filesystem>
#include <format>
#include <fstream>
#include <stdexcept>
#include <string>

namespace eonc {

std::vector<std::string> SaddleSearchJob::run() {
  if (params.saddle_search_options().method == "min_mode") {
    eonc::requirePrfoRefine(params);
  }
  std::string reactantFilename = eonc::helpers::getRelevantFile("pos.con");
  std::string displacementFilename("displacement.con");
  std::string modeFilename("direction.dat");

  if (params.main_options().checkpoint) {
    if (std::filesystem::exists("displacement_cp.con") &&
        std::filesystem::exists("mode_cp.dat")) {
      displacementFilename = "displacement_cp.con";
      modeFilename = "mode_cp.dat";
      QUILL_LOG_DEBUG(log, "Resuming from checkpoint");
    } else {
      QUILL_LOG_DEBUG(log, "No checkpoint files found");
    }
  }

  initial = std::make_shared<Matter>(pot, params);
  displacement = std::make_shared<Matter>(pot, params);
  saddle = std::make_shared<Matter>(pot, params);

  if (!eonc::io::io_ok(initial->con2matter(reactantFilename))) {
    EONC_LOG_CRITICAL("Failed to load {}", reactantFilename);
    throw std::runtime_error("failed to load " + reactantFilename);
  }

  const bool standaloneARTn = params.saddle_search_options().method == "artn";

  AtomMatrix mode = AtomMatrix::Zero(initial->numberOfAtoms(), 3);
  if (!standaloneARTn && params.saddle_search_options().displace_type ==
                             eonc::EpiCenters::DISP_LOAD) {
    // Load displacement.con, or synthesize from pos.con + direction.dat (#79).
    if (!eonc::helpers::loadOrSynthesizeDisplacement(
            *saddle, *initial, displacementFilename, modeFilename,
            params.saddle_search_options().displace_magnitude)) {
      EONC_LOG_CRITICAL("Failed to load {} (and no usable {})",
                        displacementFilename, modeFilename);
      throw std::runtime_error("missing displacement.con and direction.dat");
    }
    if (std::filesystem::exists(modeFilename)) {
      mode = eonc::helpers::loadMode(modeFilename, initial->numberOfAtoms());
    }
  } else if (!standaloneARTn && eonc::helpers::applyClientDisplacement(
                                    *saddle, *initial, params, &mode)) {
    // listed_atoms / random / last_atom / least_coordinated / not_fcc_hcp
  } else {
    *saddle = *initial;
  }

  if (standaloneARTn && std::filesystem::exists(modeFilename)) {
    mode = eonc::helpers::loadMode(modeFilename, initial->numberOfAtoms());
  }

  static_cast<void>(runPrepared(mode));
  return returnFiles;
}

std::shared_ptr<Matter>
SaddleSearchJob::runFromMatter(std::shared_ptr<Matter> seed) {
  if (!seed) {
    throw std::runtime_error("SaddleSearchJob::runFromMatter: null Matter");
  }
  initial = seed;
  initial->setPotential(pot);
  displacement = std::make_shared<Matter>(pot, params);
  saddle = std::make_shared<Matter>(pot, params);
  AtomMatrix mode = AtomMatrix::Zero(initial->numberOfAtoms(), 3);
  if (!eonc::helpers::applyClientDisplacement(*saddle, *initial, params,
                                              &mode)) {
    *saddle = *initial;
  }
  *displacement = *saddle;
  return runPrepared(mode);
}

std::shared_ptr<Matter> SaddleSearchJob::runPrepared(const AtomMatrix &mode) {
  const auto &ircOpts = params.saddle_search_options().irc;
  if (ircOpts.after) {
    if (!ircAvailable()) {
      throw std::runtime_error("[Saddle Search] irc_after needs eOn built with "
                               "-Dwith_rgsaddle=true against rgsaddle C ABI "
                               "minor 14 or newer");
    }
    if (!(ircOpts.step_size > 0.0) || !(ircOpts.converged_force > 0.0) ||
        ircOpts.max_steps < 1) {
      throw std::invalid_argument("[Saddle Search] irc_step_size and "
                                  "irc_converged_force must be positive and "
                                  "irc_max_steps at least 1");
    }
  }
  const bool useStandaloneARTn =
      params.saddle_search_options().method == "artn";
  const bool useARTnAsMinMode =
      params.saddle_search_options().method == "min_mode" &&
      params.saddle_search_options().minmode_method == "artn";

#ifdef WITH_ARTN
  if (useStandaloneARTn || useARTnAsMinMode) {
    saddleSearch =
        std::make_unique<ARTnSaddleSearch>(saddle, pot, mode, params);
  } else
#endif
  {
#ifndef WITH_ARTN
    if (useStandaloneARTn) {
      throw std::runtime_error(
          "saddle_search.method=artn requires a build with ARTn support");
    }
    if (useARTnAsMinMode) {
      throw std::runtime_error(
          "saddle_search.minmode_method=artn requires a build with ARTn "
          "support");
    }
#endif
    saddleSearch = std::make_unique<MinModeSaddleSearch>(
        saddle, mode, initial->getPotentialEnergy(), params, pot);
  }

  int status = doSaddleSearch();
  printEndState(status);
  if (ircOpts.after && status == MinModeSaddleSearch::STATUS_GOOD) {
    if (saddleSearch->getEigenvalue() < 0.0) {
      runIrc();
    } else {
      QUILL_LOG_WARNING(log,
                        "[Saddle Search] IRC skipped: lowest eigenvalue "
                        "{:f} is not negative",
                        saddleSearch->getEigenvalue());
    }
  }
  saveData(status);

  return saddle;
}

int SaddleSearchJob::doSaddleSearch() {
  long status;
  const size_t registryBefore = potRegistryForceCalls();
  const size_t potBefore = this->pot->forceCallCounter;
  try {
    status = saddleSearch->run();
  } catch (int e) {
    if (e == 100) {
      status = MinModeSaddleSearch::STATUS_POTENTIAL_FAILED;
    } else {
      QUILL_LOG_ERROR(log, "unknown exception: {}", e);
      throw e;
    }
  } catch (const std::exception &e) {
    QUILL_LOG_ERROR(log, "Saddle search potential failed: {}", e.what());
    status = MinModeSaddleSearch::STATUS_POTENTIAL_FAILED;
  }

  if (params.saddle_search_options().method == "min_mode" &&
      params.saddle_search_options().minmode_method ==
          LowestEigenmode::MINMODE_GPRDIMER) {
    fCallsSaddle = saddleSearch->getForceCalls();
  } else if (params.saddle_search_options().method == "artn") {
    fCallsSaddle = saddleSearch->getForceCalls();
  } else {
    const size_t calls = potRegistryForceCalls() - registryBefore;
    fCallsImages = calls - (this->pot->forceCallCounter - potBefore);
    fCallsSaddle += static_cast<int>(calls);
  }

  return status;
}

void SaddleSearchJob::saveData(int status) {
  std::string resultsFilename("results.dat");
  returnFiles.push_back(resultsFilename);

  std::ofstream out(resultsFilename, std::ios::binary);
  if (out) {
    out << std::format("{} termination_reason\n", status);
    out << std::format("{} termination_reason_text\n",
                       saddleSearch->describeStatus(status));
    out << "saddle_search job_type\n";
    out << std::format("{} random_seed\n", params.main_options().randomSeed);
    out << std::format(
        "{} potential_type\n",
        magic_enum::enum_name<PotType>(params.potential_options().potential));
    out << std::format("{} total_force_calls\n",
                       this->pot->forceCallCounter.load() + fCallsImages);
    out << std::format("{} force_calls_saddle\n", fCallsSaddle);
    out << std::format("{} iterations\n", saddleSearch->getIterationCount());
    if (status != MinModeSaddleSearch::STATUS_POTENTIAL_FAILED) {
      out << std::format("{:f} potential_energy_saddle\n",
                         saddle->getPotentialEnergy());
      out << std::format("{:f} final_eigenvalue\n",
                         saddleSearch->getEigenvalue());
    }
    out << std::format("{:f} potential_energy_reactant\n",
                       initial->getPotentialEnergy());
    if (irc.ran) {
      out << std::format("{} irc_branches_at_minimum\n", irc.branchesAtMinimum);
      out << std::format("{} irc_ends_relaxed\n", int{irc.endsRelaxed});
      out << std::format("{} irc_connects_reactant\n",
                         int{irc.connectsReactant});
      out << std::format("{} force_calls_irc\n", irc.forceCalls);
      out << std::format("{} force_calls_irc_minimization\n",
                         irc.forceCallsMin);
      out << std::format("{:f} potential_energy_irc_min1\n", irc.energies[0]);
      out << std::format("{:f} potential_energy_irc_min2\n", irc.energies[1]);
    }
  }

  std::string modeFilename("mode.dat");
  returnFiles.push_back(modeFilename);
  eonc::helpers::saveMode(modeFilename, saddle, saddleSearch->getEigenvector());

  std::string saddleFilename("saddle.con");
  returnFiles.push_back(saddleFilename);
  if (!eonc::io::io_ok(saddle->matter2con(saddleFilename))) {
    QUILL_LOG_ERROR(log, "Failed to write {}", saddleFilename);
  }
}

void SaddleSearchJob::runIrc() {
  irc.ran = true;
  const size_t before = pot->forceCallCounter;
  IrcPath path;
  try {
    path = followIrc(*saddle, saddleSearch->getEigenvector(), params);
  } catch (const std::exception &e) {
    irc.forceCalls = pot->forceCallCounter - before;
    QUILL_LOG_ERROR(log, "[Saddle Search] IRC failed: {}", e.what());
    return;
  }
  irc.forceCalls = pot->forceCallCounter - before;
  irc.branchesAtMinimum =
      int{path.minus.reachedMinimum} + int{path.plus.reachedMinimum};

  // irc.con runs from the minus end through the saddle to the plus end.
  // irc.dat follows the neb.dat layout, with max_force in the last column.
  // rxn_coord and energy are taken from the saddle.
  struct Point {
    const AtomMatrix *pos;
    double energy;
    double maxForce;
  };
  std::vector<Point> points;
  for (size_t i = path.minus.positions.size(); i-- > 0;) {
    points.push_back({&path.minus.positions[i], path.minus.energies[i],
                      path.minus.maxForces[i]});
  }
  const size_t saddleIndex = points.size();
  const AtomMatrix saddlePos = saddle->getPositions();
  const double saddleEnergy = saddle->getPotentialEnergy();
  points.push_back({&saddlePos, saddleEnergy, saddle->maxForce()});
  for (size_t i = 0; i < path.plus.positions.size(); ++i) {
    points.push_back({&path.plus.positions[i], path.plus.energies[i],
                      path.plus.maxForces[i]});
  }
  std::vector<double> arc(points.size(), 0.0);
  for (size_t i = 1; i < points.size(); ++i) {
    arc[i] =
        arc[i - 1] + saddle->pbc(*points[i].pos - *points[i - 1].pos).norm();
  }
  Matter frame(*saddle);
  std::ofstream dat("irc.dat", std::ios::binary);
  dat << std::format("{:>3s} {:>12s} {:>12s} {:>12s}\n", "img", "rxn_coord",
                     "energy", "max_force");
  for (size_t i = 0; i < points.size(); ++i) {
    frame.setPositions(*points[i].pos);
    eonc::io::ConFrameMetadata meta;
    meta.energy = points[i].energy;
    if (!eonc::io::io_ok(frame.matter2con("irc.con", i > 0, &meta))) {
      QUILL_LOG_ERROR(log, "Failed to write irc.con");
    }
    dat << std::format("{:>3} {:>12.6f} {:>12.6f} {:>12.6f}\n", i,
                       arc[i] - arc[saddleIndex],
                       points[i].energy - saddleEnergy, points[i].maxForce);
  }
  returnFiles.push_back("irc.con");
  returnFiles.push_back("irc.dat");

  const size_t minBefore = pot->forceCallCounter;
  irc.endsRelaxed = true;
  const IrcBranch *branches[2] = {&path.minus, &path.plus};
  for (int k = 0; k < 2; ++k) {
    Matter end(*saddle);
    if (!branches[k]->positions.empty()) {
      end.setPositions(branches[k]->positions.back());
    }
    irc.endsRelaxed = end.relax(true) && irc.endsRelaxed;
    irc.energies[k] = end.getPotentialEnergy();
    Matter reactant(*initial);
    irc.connectsReactant = irc.connectsReactant || reactant.compare(end);
    const std::string name = std::format("irc_min{}.con", k + 1);
    returnFiles.push_back(name);
    if (!eonc::io::io_ok(end.matter2con(name))) {
      QUILL_LOG_ERROR(log, "Failed to write {}", name);
    }
  }
  irc.forceCallsMin = pot->forceCallCounter - minBefore;
  QUILL_LOG_INFO(log,
                 "[Saddle Search] IRC: {} of 2 branches at a minimum, {} force "
                 "calls, {} more to relax the ends",
                 irc.branchesAtMinimum, irc.forceCalls, irc.forceCallsMin);
}

void SaddleSearchJob::printEndState(int status) {
  auto msg = saddleSearch->describeStatus(status);
  if (status == MinModeSaddleSearch::STATUS_GOOD) {
    QUILL_LOG_DEBUG(log, "[Saddle Search] {}", msg);
  } else {
    QUILL_LOG_WARNING(log, "[Saddle Search] {}", msg);
  }
}

} // namespace eonc
