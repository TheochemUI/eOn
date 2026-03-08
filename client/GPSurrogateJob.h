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
#include "PyGuard.h"
#include <format>
#include <pybind11/eigen.h>

#include "NudgedElasticBand.h"

namespace eonc {


class GPSurrogateJob : public Job {
public:
  GPSurrogateJob(std::unique_ptr<Parameters> parameters)
      : Job(std::move(parameters)) {
    eonc::ensure_interpreter();
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
};

namespace helpers::surrogate {
MatrixXd get_features(const std::vector<Matter> &matobjs);
MatrixXd get_features(const std::vector<std::shared_ptr<Matter>> &matobjs);
MatrixXd get_targets(std::vector<std::shared_ptr<Matter>> &matobjs,
                     std::shared_ptr<Potential> true_pot);
MatrixXd get_targets(std::vector<Matter> &matobjs,
                     std::shared_ptr<Potential> true_pot);
Eigen::VectorXd make_target(Matter &m1, std::shared_ptr<Potential> true_pot);
std::pair<Eigen::VectorXd, Eigen::VectorXd>
getNewDataPoint(const std::vector<std::shared_ptr<Matter>> &matobjs,
                std::shared_ptr<Potential> true_pot);
std::vector<Matter> getMidSlice(const std::vector<Matter> &matobjs);
bool accuratePES(std::vector<std::shared_ptr<Matter>> &matobjs,
                 std::shared_ptr<Potential> true_pot);
std::pair<double, Eigen::VectorXd::Index>
getMaxUncertainty(const std::vector<std::shared_ptr<Matter>> &matobjs);
} // namespace helpers::surrogate

namespace helpers::eigen {
MatrixXd vertCat(const MatrixXd &m1, const MatrixXd &m2);
void addVectorRow(MatrixXd &data, const Eigen::VectorXd &newrow);
// Modifies data
} // namespace helpers::eigen

} // namespace eonc

using eonc::GPSurrogateJob;
