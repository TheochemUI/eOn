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
#include "HelperFunctions.h"
#include "Hessian.h"
#include "Matter.h"
#include "Potential.h"
#include "Prefactor.h"

#include <cmath>
#include <format>
#include <fstream>
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

  reactant->con2matter("reactant.con");
  saddle->con2matter("saddle.con");
  product->con2matter("product.con");
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
    reactant->con2matter(matterFilename);
    saddle->con2matter(matterFilename);
    product->con2matter(matterFilename);

    atoms = eonc::Prefactor::allFreeAtoms(reactant.get());
  } else {
    reactant->con2matter(reactantFilename);
    saddle->con2matter(saddleFilename);
    product->con2matter(productFilename);

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
