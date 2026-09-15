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
#include "eon/PrefactorJob.h"
#include "eon/BaseStructures.h"
#include "eon/EonLogger.h"
#include "eon/HelperFunctions.h"
#include "eon/Hessian.h"
#include "eon/JobResult.h"
#include "eon/Matter.h"
#include "eon/PotRegistry.h"
#include "eon/Potential.h"
#include "eon/Prefactor.h"

#include <cmath>
#include <filesystem>
#include <format>
#include <fstream>
#include <stdexcept>
#include <string>

namespace eonc {

const char PrefactorJob::PREFACTOR_REACTANT[] = "reactant";
const char PrefactorJob::PREFACTOR_SADDLE[] = "saddle";
const char PrefactorJob::PREFACTOR_PRODUCT[] = "product";

std::vector<std::string> PrefactorJob::run() {
  std::vector<std::string> returnFiles;
  VectorXd freqs;

  std::string reactantFilename = eonc::helpers::getRelevantFile("reactant.con");
  std::string saddleFilename = eonc::helpers::getRelevantFile("saddle.con");
  std::string productFilename = eonc::helpers::getRelevantFile("product.con");

  auto reactant = std::make_unique<Matter>(pot, params);
  auto saddle = std::make_unique<Matter>(pot, params);
  auto product = std::make_unique<Matter>(pot, params);

  if (!eonc::io::io_ok(reactant->con2matter(reactantFilename)) ||
      !eonc::io::io_ok(saddle->con2matter(saddleFilename)) ||
      !eonc::io::io_ok(product->con2matter(productFilename))) {
    EONC_LOG_CRITICAL("Failed to load reactant/saddle/product for prefactor");
    throw std::runtime_error("failed to load prefactor geometries");
  }
  double pref1 = 0.0, pref2 = 0.0;
  const int prefStatus = eonc::Prefactor::getPrefactors(
      params, reactant.get(), saddle.get(), product.get(), pref1, pref2);

  VectorXi atoms;
  if (params.prefactor_options.all_free_atoms) {
    std::string matterFilename;
    if (params.prefactor_options.configuration ==
        PrefactorJob::PREFACTOR_REACTANT) {
      matterFilename = reactantFilename;
    } else if (params.prefactor_options.configuration ==
               PrefactorJob::PREFACTOR_SADDLE) {
      matterFilename = saddleFilename;
    } else if (params.prefactor_options.configuration ==
               PrefactorJob::PREFACTOR_PRODUCT) {
      matterFilename = productFilename;
    }
    if (!eonc::io::io_ok(reactant->con2matter(matterFilename)) ||
        !eonc::io::io_ok(saddle->con2matter(matterFilename)) ||
        !eonc::io::io_ok(product->con2matter(matterFilename))) {
      EONC_LOG_CRITICAL("Failed to reload {} for all-free-atoms prefactor",
                        matterFilename);
      throw std::runtime_error("failed to load prefactor configuration");
    }

    atoms = eonc::Prefactor::allFreeAtoms(reactant.get());
  } else {
    if (!eonc::io::io_ok(reactant->con2matter(reactantFilename)) ||
        !eonc::io::io_ok(saddle->con2matter(saddleFilename)) ||
        !eonc::io::io_ok(product->con2matter(productFilename))) {
      EONC_LOG_CRITICAL(
          "Failed to reload reactant/saddle/product for prefactor");
      throw std::runtime_error("failed to load prefactor geometries");
    }

    atoms = eonc::Prefactor::movedAtoms(params, reactant.get(), saddle.get(),
                                        product.get());
  }
  bool failed = (prefStatus == -1) || (atoms.rows() == 0);

  if (!failed) {
    if (params.prefactor_options.configuration ==
        PrefactorJob::PREFACTOR_REACTANT) {
      Hessian hessian(params, reactant.get());
      freqs = hessian.getFreqs(reactant.get(), atoms);
    } else if (params.prefactor_options.configuration ==
               PrefactorJob::PREFACTOR_SADDLE) {
      Hessian hessian(params, saddle.get());
      freqs = hessian.getFreqs(saddle.get(), atoms);
    } else if (params.prefactor_options.configuration ==
               PrefactorJob::PREFACTOR_PRODUCT) {
      Hessian hessian(params, product.get());
      freqs = hessian.getFreqs(product.get(), atoms);
    }
  }

  if (!failed) {
    failed = freqs.size() != 3 * atoms.rows();
  }

  std::string results_file("results.dat");
  std::string freq_file("freq.dat");
  returnFiles.push_back(results_file);
  returnFiles.push_back(freq_file);

  std::ofstream outFreq(freq_file, std::ios::binary);

  auto env = JobResultEnvelope::fromMinimization(
      failed ? RunStatus::FAIL_POTENTIAL_FAILED : RunStatus::GOOD,
      params.potential_options.potential,
      PotRegistry::get().total_force_calls(), false, 0.0);
  env.job_type = "prefactor";
  env.tags.emplace_back("good", failed ? "false" : "true");
  env.extras.emplace_back("force_calls", static_cast<double>(env.force_calls));
  if (!failed) {
    env.extras.emplace_back("prefactor_reactant_to_product", pref1);
    env.extras.emplace_back("prefactor_product_to_reactant", pref2);
  }
  env.writeResultsDat(results_file);

  if (outFreq && !failed) {
    for (int i = 0; i < freqs.size(); i++) {
      if (0. < freqs[i]) {
        outFreq << std::format("{:f}\n",
                               std::sqrt(freqs[i]) /
                                   (2 * eonc::helpers::pi * 10.18e-15));
      } else {
        outFreq << std::format("{:f}\n",
                               -std::sqrt(-freqs[i]) /
                                   (2 * eonc::helpers::pi * 10.18e-15));
      }
    }
  }

  if (std::filesystem::exists("freqs.dat")) {
    returnFiles.push_back("freqs.dat");
  }
  if (std::filesystem::exists("hessian.dat")) {
    returnFiles.push_back("hessian.dat");
  }

  return returnFiles;
}

} // namespace eonc
