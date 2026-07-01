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
#include "PrefactorJob.h"
#include "EonLogger.h"
#include "HelperFunctions.h"
#include "Hessian.h"
#include "Matter.h"
#include "Potential.h"
#include "Prefactor.h"

#include <cmath>
#include <format>
#include <fstream>
#include <stdexcept>
#include <string>

const char PrefactorJob::PREFACTOR_REACTANT[] = "reactant";
const char PrefactorJob::PREFACTOR_SADDLE[] = "saddle";
const char PrefactorJob::PREFACTOR_PRODUCT[] = "product";

std::vector<std::string> PrefactorJob::run() {
  std::vector<std::string> returnFiles;
  VectorXd freqs;

  std::string reactantFilename("reactant.con");
  std::string saddleFilename("saddle.con");
  std::string productFilename("product.con");

  auto reactant = std::make_unique<Matter>(pot, params);
  auto saddle = std::make_unique<Matter>(pot, params);
  auto product = std::make_unique<Matter>(pot, params);

  if (!eonc::io::io_ok(reactant->con2matter("reactant.con")) ||
      !eonc::io::io_ok(saddle->con2matter("saddle.con")) ||
      !eonc::io::io_ok(product->con2matter("product.con"))) {
    EONC_LOG_CRITICAL("Failed to load reactant/saddle/product for prefactor");
    throw std::runtime_error("failed to load prefactor geometries");
  }
  double pref1, pref2;
  eonc::Prefactor::getPrefactors(params, reactant.get(), saddle.get(),
                                 product.get(), pref1, pref2);

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
  assert(3 * atoms.rows() > 0);

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

  bool failed = freqs.size() != 3 * atoms.rows();

  std::string results_file("results.dat");
  std::string freq_file("freq.dat");
  returnFiles.push_back(results_file);
  returnFiles.push_back(freq_file);

  std::ofstream outResults(results_file, std::ios::binary);
  std::ofstream outFreq(freq_file, std::ios::binary);

  if (outResults) {
    outResults << std::format("{} good\n", failed ? "false" : "true");
    outResults << std::format("{} force_calls\n",
                              PotRegistry::get().total_force_calls());
  }

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

  return returnFiles;
}
