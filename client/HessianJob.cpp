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
#include "HessianJob.h"
#include "Hessian.h"
#include "Matter.h"
#include "Potential.h"

#include <format>
#include <fstream>
#include <string>

std::vector<std::string> HessianJob::run(void) {
  std::string matter_in("pos.con");

  std::vector<std::string> returnFiles;

  auto matter = std::make_unique<Matter>(pot, params);

  matter->con2matter(matter_in);

  Hessian hessian(params, matter.get());
  long nAtoms = matter->numberOfAtoms();

  VectorXi moved(nAtoms);
  moved.setConstant(-1);

  int nMoved = 0;
  // [Hessian] atom_list = mobile/displaced atoms for FD (hybrid/PHVA-class
  // active set), comma-separated 0-based indices, or "All" = every non-fixed
  // atom. Intersect with free flags so frozen CON atoms are never displaced.
  const std::string &atomList = params.hessian_options.atom_list;
  const bool useExplicitList =
      !atomList.empty() &&
      atomList != "All" && atomList != "all" && atomList != "ALL";

  if (useExplicitList) {
    std::string token;
    for (size_t p = 0; p <= atomList.size(); ++p) {
      const char c = (p < atomList.size()) ? atomList[p] : ',';
      if (c == ',' || c == ' ' || c == '\t' || p == atomList.size()) {
        if (!token.empty()) {
          try {
            const long idx = std::stol(token);
            if (idx >= 0 && idx < nAtoms && !matter->getFixed(idx)) {
              moved[nMoved++] = static_cast<int>(idx);
            }
          } catch (const std::exception &) {
            // skip non-integer tokens
          }
          token.clear();
        }
      } else {
        token.push_back(c);
      }
    }
  } else {
    for (int i = 0; i < nAtoms; i++) {
      if (!matter->getFixed(i)) {
        moved[nMoved] = i;
        nMoved++;
      }
    }
  }
  moved = moved.head(nMoved);
  if (nMoved == 0) {
    // No free atoms: leave results.dat with force_calls only (no crash).
    std::string results_file("results.dat");
    returnFiles.push_back(results_file);
    std::ofstream out(results_file, std::ios::binary);
    if (out) {
      out << std::format("{} force_calls\n",
                         PotRegistry::get().total_force_calls());
    }
    return returnFiles;
  }
  hessian.getFreqs(matter.get(), moved);

  std::string results_file("results.dat");
  returnFiles.push_back(results_file);

  std::ofstream out(results_file, std::ios::binary);
  if (out) {
    out << std::format("{} force_calls\n",
                       PotRegistry::get().total_force_calls());
  }

  return returnFiles;
}
