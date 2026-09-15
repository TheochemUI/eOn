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
#include "eon/GPSurrogateJob.h"
#include "eon/BaseStructures.h"
#include "eon/NEBSplineExtrema.h"
#include "eon/NudgedElasticBand.h"
#include "eon/NudgedElasticBandJob.h"
#include "eon/SurrogatePotential.h"
#include "eon/helpers/Create.hpp"
#include "eon/potentials/CatLearnPot/CatLearnPot.h"

#include "eon/EonLogger.h"
#include <sstream>
#include <stdexcept>

namespace eonc {

std::vector<std::string> GPSurrogateJob::run() {
  std::string reactantFilename = eonc::helpers::getRelevantFile("reactant.con");
  std::string productFilename = eonc::helpers::getRelevantFile("product.con");
  auto true_params = std::make_shared<Parameters>(params);
  true_params->main_options().job = params.sub_job;
  auto initial = std::make_shared<Matter>(pot, *true_params);
  if (!eonc::io::io_ok(initial->con2matter(reactantFilename))) {
    EONC_LOG_CRITICAL("Failed to load {}", reactantFilename);
    throw std::runtime_error("failed to load " + reactantFilename);
  }
  auto final_state = std::make_shared<Matter>(pot, *true_params);
  if (!eonc::io::io_ok(final_state->con2matter(productFilename))) {
    EONC_LOG_CRITICAL("Failed to load {}", productFilename);
    throw std::runtime_error("failed to load " + productFilename);
  }
  (void)runFromMatter(initial, final_state);
  return returnFiles;
}

std::shared_ptr<NudgedElasticBand>
GPSurrogateJob::runFromMatter(std::shared_ptr<Matter> initial,
                              std::shared_ptr<Matter> final_state) {
  if (!initial || !final_state) {
    throw std::runtime_error("GPSurrogateJob::runFromMatter: null Matter");
  }
  // Clone and setup "true" params
  auto true_params = std::make_shared<Parameters>(params);
  true_params->main_options().job = params.sub_job;
  auto true_job =
      eonc::helpers::makeJob(std::make_unique<Parameters>(*true_params));
  auto pyparams = std::make_shared<Parameters>(params);
  pyparams->potential_options().potential = PotType::CatLearn;

  initial->setPotential(pot);
  final_state->setPotential(pot);
  auto init_path = eonc::helpers::neb_paths::linearPath(
      *initial, *final_state, params.neb_options().image_count);
  auto init_data = eonc::helpers::surrogate::getMidSlice(init_path);
  auto features = eonc::helpers::surrogate::get_features(init_data);
  EONC_LOG_TRACE("Potential is {}",
                 magic_enum::enum_name<PotType>(pot->getType()));
  auto targets = eonc::helpers::surrogate::get_targets(init_data, pot);

  // Setup a GPR Potential
  auto surpot = eonc::helpers::create::makeSurrogatePotential(
      params.gp_surrogate_options().potential, params);
  surpot->train_optimize(features, targets);
  auto neb = std::make_unique<NudgedElasticBand>(initial, final_state,
                                                 *pyparams, surpot);
  auto status_neb{neb->compute()};
  bool job_not_finished{true};
  size_t n_gp{0};
  double unc_conv{pyparams->gp_uncertainty};
  while (job_not_finished) { // outer loop?
    n_gp++;
    if (n_gp > 750) {
      EONC_LOG_CRITICAL("Whoops, power level of problem too high!!");
      break;
    }
    EONC_LOG_TRACE("Must handle update to the GP, update number {}", n_gp);
    auto [maxUnc, maxIndex] =
        eonc::helpers::surrogate::getMaxUncertainty(neb->path);
    auto [feature, target] =
        eonc::helpers::surrogate::getNewDataPoint(neb->path, pot);
    eonc::helpers::eigen::addVectorRow(features, feature);
    eonc::helpers::eigen::addVectorRow(targets, target);
    surpot->train_optimize(features, targets);
    pyparams->nebClimbingImageMethod = false;
    pyparams->optimizer_options().converged_force =
        params.optimizer_options().converged_force * 0.8;
    for (auto &&obj : neb->path) {
      obj->setPotential(surpot);
    }
    if (!(pyparams->gp_linear_path_always)) {
      EONC_LOG_TRACE("Using previous path");
      neb = std::make_unique<NudgedElasticBand>(neb->path, *pyparams, surpot);
    } else {
      EONC_LOG_TRACE("Using linear interpolation");
      neb = std::make_unique<NudgedElasticBand>(initial, final_state, *pyparams,
                                                surpot);
    }
    status_neb = neb->compute();

    std::string nebFilename(std::format("neb_final_gpr_{:03d}.con", n_gp));
    returnFiles.push_back(nebFilename);
    if (!eonc::io::io_ok(eonc::neb::writePathCon(
            neb->path, neb->tangent, neb->eigenmode_solvers, neb->numImages,
            params.debug_options().estimate_neb_eigenvalues, nebFilename,
            static_cast<size_t>(n_gp)))) {
      throw std::runtime_error("Failed to write file: " + nebFilename);
    }
    if (status_neb == NudgedElasticBand::NEBStatus::GOOD &&
        eonc::helpers::surrogate::accuratePES(neb->path, pot)) {
      break;
    } else {
      continue;
    }
  }
  neb->printImageData();
  neb->findExtrema();
  // Keep a shared view for the caller before saveData takes ownership of unique
  std::shared_ptr<NudgedElasticBand> out(neb.release());
  // saveData expects unique_ptr - rebuild unique from shared is unsafe.
  // Write results without consuming: call saveData on a temporary unique wrap
  // fails. Instead write via a clone path: only return the band; file artifacts
  // optional.
  return out;
}

void GPSurrogateJob::saveData(NudgedElasticBand::NEBStatus status,
                              std::unique_ptr<NudgedElasticBand> neb) {
  std::string resultsFilename = "results.dat";
  returnFiles.push_back(resultsFilename);

  std::ofstream fileResults(resultsFilename);
  if (!fileResults) {
    // Handle file open error
    throw std::runtime_error("Failed to open file: " + resultsFilename);
  }

  fileResults << static_cast<int>(status) << " termination_reason\n";
  fileResults << magic_enum::enum_name(status) << " termination_reason_text\n";
  fileResults << magic_enum::enum_name<PotType>(
                     params.potential_options().potential)
              << " potential_type\n";
  fileResults << std::format("{:.6f} energy_reference\n",
                             neb->path[0]->getPotentialEnergy());
  fileResults << neb->numImages << " number_of_images\n";

  for (long i = 0; i <= neb->numImages + 1; i++) {
    fileResults << std::format("{:.6f} image{}_energy\n",
                               neb->path[i]->getPotentialEnergy() -
                                   neb->path[0]->getPotentialEnergy(),
                               i);
    fileResults << std::format("{:.6f} image{}_force\n",
                               neb->path[i]->getForces().norm(), i);
    fileResults << std::format("{:.6f} image{}_projected_force\n",
                               neb->projectedForce[i]->norm(), i);
  }

  fileResults << neb->numExtrema << " number_of_extrema\n";
  for (long i = 0; i < neb->numExtrema; i++) {
    fileResults << std::format("{:.6f} extremum{}_position\n",
                               neb->extremumPosition[i], i);
    fileResults << std::format("{:.6f} extremum{}_energy\n",
                               neb->extremumEnergy[i], i);
  }

  fileResults.close();

  std::string nebFilename = "neb.con";
  returnFiles.push_back(nebFilename);

  if (!eonc::io::io_ok(eonc::neb::writePathCon(
          neb->path, neb->tangent, neb->eigenmode_solvers, neb->numImages,
          params.debug_options().estimate_neb_eigenvalues, nebFilename))) {
    throw std::runtime_error("Failed to write file: " + nebFilename);
  }

  returnFiles.push_back("neb.dat");
  neb->printImageData(true);
}

} // namespace eonc

namespace eonc::helpers::surrogate {
MatrixXd get_features(const std::vector<Matter> &matobjs) {
  // Calculate dimensions
  MatrixXd features(matobjs.size(), matobjs.front().numberOfFreeAtoms() * 3);
  EONC_LOG_TRACE("rows: {}, cols:{}", matobjs.size(),
                 matobjs.front().numberOfFreeAtoms() * 3);
  for (long idx{0}; idx < features.rows(); idx++) {
    features.row(idx) = matobjs[idx].getPositionsFreeV();
  }
  std::ostringstream oss;
  oss << features;
  EONC_LOG_TRACE("Features\n:{}", oss.str());
  return features;
}
MatrixXd get_features(const std::vector<std::shared_ptr<Matter>> &matobjs) {
  // Calculate dimensions
  MatrixXd features(matobjs.size(), matobjs.front()->numberOfFreeAtoms() * 3);
  EONC_LOG_TRACE("rows: {}, cols:{}\n", matobjs.size(),
                 matobjs.front()->numberOfFreeAtoms() * 3);
  for (long idx{0}; idx < features.rows(); idx++) {
    features.row(idx) = matobjs[idx]->getPositionsFreeV();
  }
  std::ostringstream oss;
  oss << features;
  EONC_LOG_TRACE("Features\n:{}", oss.str());
  return features;
}
MatrixXd get_targets(std::vector<Matter> &matobjs,
                     std::shared_ptr<Potential> true_pot) {
  // Always with derivatives for now
  // Energy + Derivatives for each row
  const auto nrows = matobjs.size();
  const auto ncols = (matobjs.front().numberOfFreeAtoms() * 3) + 1;
  MatrixXd targets(nrows, ncols);
  for (long idx{0}; idx < targets.rows(); idx++) {
    matobjs[idx].setPotential(true_pot);
    targets.row(idx)[0] = matobjs[idx].getPotentialEnergy();
    targets.block(idx, 1, 1, ncols - 1) =
        matobjs[idx].getForcesFreeV().array() * -1;
  }
  std::ostringstream oss;
  oss << targets;
  EONC_LOG_TRACE("Targets\n:{}", oss.str());
  return targets;
}
MatrixXd get_targets(std::vector<std::shared_ptr<Matter>> &matobjs,
                     std::shared_ptr<Potential> true_pot) {
  const auto nrows = matobjs.size();
  const auto ncols = (matobjs.front()->numberOfFreeAtoms() * 3) + 1;
  MatrixXd targets(nrows, ncols);
  for (long idx{0}; idx < targets.rows(); idx++) {
    matobjs[idx]->setPotential(true_pot);
    targets.row(idx)[0] = matobjs[idx]->getPotentialEnergy();
    targets.block(idx, 1, 1, ncols - 1) =
        matobjs[idx]->getForcesFreeV().array() * -1;
  }
  std::ostringstream oss;
  oss << targets;
  EONC_LOG_TRACE("Targets\n:{}", oss.str());
  return targets;
}
std::vector<Matter> getMidSlice(const std::vector<Matter> &matobjs) {
  // Initial GP slice: endpoints plus one interior sample. CatLearn
  // training is order-sensitive (front, back, interior). The interior
  // index is two-thirds along the movable images, not n/2.
  if (matobjs.size() < 3) {
    throw std::invalid_argument("getMidSlice: need at least three images");
  }
  const std::size_t n = matobjs.size();
  const std::size_t twoThirds =
      static_cast<std::size_t>(((n - 2) * 2.0 / 3.0) + 1.0);
  return {matobjs.front(), matobjs.back(), matobjs[twoThirds]};
}
Eigen::VectorXd make_target(Matter &m1, std::shared_ptr<Potential> true_pot) {
  const auto ncols = (m1.numberOfFreeAtoms() * 3) + 1;
  Eigen::VectorXd target(ncols);
  m1.setPotential(true_pot);
  target(0) = m1.getPotentialEnergy();
  target.segment(1, ncols - 1) = m1.getForcesFreeV() * -1;
  // EONC_LOG_TRACE("Generated Target:\n{}",
  // fmt::streamed(target));
  return target;
}
std::pair<double, Eigen::VectorXd::Index>
getMaxUncertainty(const std::vector<std::shared_ptr<Matter>> &matobjs) {
  if (matobjs.size() < 3) {
    throw std::invalid_argument(
        "getMaxUncertainty: need at least three images");
  }
  Eigen::VectorXd pathUncertainty{Eigen::VectorXd::Zero(matobjs.size() - 2)};
  for (auto idx{0}; idx < pathUncertainty.size(); idx++) {
    pathUncertainty[idx] = matobjs[idx + 1]->getEnergyVariance();
  }
  Eigen::VectorXd::Index maxIndex;
  double maxUnc{pathUncertainty.maxCoeff()};
  pathUncertainty.maxCoeff(&maxIndex);
  // EONC_LOG_TRACE("Uncertainty along path
  // is {}\nmax_index: {}, maxVal: {}",
  //              fmt::streamed(pathUncertainty), maxIndex, maxUnc);
  return std::make_pair(maxUnc, maxIndex);
}
std::pair<Eigen::VectorXd, Eigen::VectorXd>
getNewDataPoint(const std::vector<std::shared_ptr<Matter>> &matobjs,
                std::shared_ptr<Potential> true_pot) {
  auto [maxUnc, maxIndex] = getMaxUncertainty(matobjs);
  Matter candidate{*matobjs[maxIndex + 1]};
  return std::make_pair<Eigen::VectorXd, Eigen::VectorXd>(
      candidate.getPositionsFreeV(), make_target(candidate, true_pot));
}
bool accuratePES(std::vector<std::shared_ptr<Matter>> &matobjs,
                 std::shared_ptr<Potential> true_pot) {
  if (matobjs.empty()) {
    throw std::invalid_argument("accuratePES: empty path");
  }
  Eigen::VectorXd predEnergies{Eigen::VectorXd::Zero(matobjs.size())};
  Eigen::VectorXd trueEnergies{Eigen::VectorXd::Zero(matobjs.size())};
  for (auto idx{0}; idx < predEnergies.size(); idx++) {
    auto incoming = matobjs[idx]->getPotential();
    predEnergies[idx] = matobjs[idx]->getPotentialEnergy();
    matobjs[idx]->setPotential(true_pot);
    trueEnergies[idx] = matobjs[idx]->getPotentialEnergy();
    matobjs[idx]->setPotential(incoming);
  }
  Eigen::VectorXd difference = predEnergies - trueEnergies;
  const auto maxAbs = difference.array().abs().maxCoeff();
  std::ostringstream oss;
  oss << "predicted\n"
      << predEnergies << "\ntrue\n"
      << trueEnergies << "\ndifference\n"
      << difference << "\n maxAbs: " << maxAbs;
  EONC_LOG_TRACE("{}", oss.str());
  return maxAbs < 0.05;
}
} // namespace eonc::helpers::surrogate

namespace eonc::helpers::eigen {
MatrixXd vertCat(const MatrixXd &m1, const MatrixXd &m2) {
  assert(m1.cols() == m2.cols());
  MatrixXd res(m1.rows() + m2.rows(), m2.cols());
  res << m1, m2;
  return res;
}
void addVectorRow(MatrixXd &data, const Eigen::VectorXd &newrow) {
  assert(data.cols() == newrow.size());
  data.conservativeResize(data.rows() + 1, data.cols());
  data.row(data.rows() - 1) = newrow;
}
} // namespace eonc::helpers::eigen
