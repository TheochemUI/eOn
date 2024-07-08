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
#pragma once

#include "HelperFunctions.h"
#include "Job.h"
#include "Parameters.h"

#ifdef WITH_CATLEARN
#include "potentials/CatLearnPot/CatLearnPot.h"
#endif
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <pybind11/eigen.h>
#include <pybind11/embed.h>

#include "NudgedElasticBand.h"

class GPSurrogateJob : public Job {
public:
  GPSurrogateJob(std::unique_ptr<Parameters> parameters)
      : Job(std::move(parameters)) {
    // debugging
#ifndef NDEBUG
    py::module_ sys_mod = py::module_::import("sys");
    py::module_ ipdb_mod = py::module_::import("ipdb");
    sys_mod.attr("breakpointhook") = ipdb_mod.attr("set_trace");
#endif // NDEBUG
  }
  ~GPSurrogateJob(void) = default;
  std::vector<std::string> run(void) override;

private:
  void saveData(NudgedElasticBand::NEBStatus status,
                std::unique_ptr<NudgedElasticBand> neb);
  std::vector<std::string> returnFiles;
#ifndef WITH_ASE_ORCA
  // When ASE_ORCA is used as a potential ClientEON.cpp holds lock in main
  pybind11::scoped_interpreter guard{};
#endif
};

namespace helper_functions::surrogate {
MatrixType get_features(const std::vector<Matter> &matobjs);
MatrixType get_features(const std::vector<std::shared_ptr<Matter>> &matobjs);
MatrixType get_targets(std::vector<std::shared_ptr<Matter>> &matobjs,
                       std::shared_ptr<Potential> true_pot);
MatrixType get_targets(std::vector<Matter> &matobjs,
                       std::shared_ptr<Potential> true_pot);
VectorType make_target(Matter &m1, std::shared_ptr<Potential> true_pot);
std::pair<VectorType, VectorType>
getNewDataPoint(const std::vector<std::shared_ptr<Matter>> &matobjs,
                std::shared_ptr<Potential> true_pot);
std::vector<Matter> getMidSlice(const std::vector<Matter> &matobjs);
bool accuratePES(std::vector<std::shared_ptr<Matter>> &matobjs,
                 std::shared_ptr<Potential> true_pot);
std::pair<double, VectorType::Index>
getMaxUncertainty(const std::vector<std::shared_ptr<Matter>> &matobjs);
} // namespace helper_functions::surrogate

namespace helper_functions::eigen {
MatrixType vertCat(const MatrixType &m1, const MatrixType &m2);
void addVectorRow(MatrixType &data, const VectorType &newrow);
// Modifies data
} // namespace helper_functions::eigen
