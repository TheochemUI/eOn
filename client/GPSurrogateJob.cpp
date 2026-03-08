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
#include "GPSurrogateJob.h"
#include "BaseStructures.h"
#include "NudgedElasticBand.h"
#include "NudgedElasticBandJob.h"
#include "SurrogatePotential.h"
#include "helpers/Create.hpp"
#include "potentials/CatLearnPot/CatLearnPot.h"

#include "EonLogger.h"
#include <sstream>
using namespace std;
std::vector<std::string> GPSurrogateJob::run(void) {
  // Start working
  std::string reactantFilename =
      eonc::helpers::getRelevantFile("reactant.con");
  std::string productFilename =
      eonc::helpers::getRelevantFile("product.con");

  // Clone and setup "true" params
  auto true_params = std::make_shared<Parameters>(params);
  true_params->main_options.job = params.sub_job;
  auto true_job =
      eonc::helpers::makeJob(std::make_unique<Parameters>(*true_params));
  auto pyparams = std::make_shared<Parameters>(params);
  pyparams->potential_options.potential = PotType::CatLearn;

  // Get possible initial data source
  auto initial = std::make_shared<Matter>(pot, *true_params);
  initial->con2matter(reactantFilename);
  auto final_state = std::make_shared<Matter>(pot, *true_params);
  final_state->con2matter(productFilename);
  auto init_path = eonc::helpers::neb_paths::linearPath(
      *initial, *final_state, params.neb_options.image_count);
  auto init_data = eonc::helpers::surrogate::getMidSlice(init_path);
  auto features = eonc::helpers::surrogate::get_features(init_data);
  EONC_LOG_TRACE("Potential is {}",
                 magic_enum::enum_name<PotType>(pot->getType()));
  auto targets = eonc::helpers::surrogate::get_targets(init_data, pot);

  // Setup a GPR Potential
  auto surpot = eonc::helpers::create::makeSurrogatePotential(
      params.gp_surrogate_options.potential, params);
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
    // if (status_neb == NudgedElasticBand::NEBStatus::MAX_UNCERTAINTY ||
    // status_neb == NudgedElasticBand::NEBStatus::BAD_MAX_ITERATIONS) {
    EONC_LOG_TRACE("Must handle update to the GP, update number {}", n_gp);
    auto [maxUnc, maxIndex] =
        eonc::helpers::surrogate::getMaxUncertainty(neb->path);
    auto [feature, target] =
        eonc::helpers::surrogate::getNewDataPoint(neb->path, pot);
    eonc::helpers::eigen::addVectorRow(features, feature);
    eonc::helpers::eigen::addVectorRow(targets, target);
    surpot->train_optimize(features, targets);
    pyparams->nebClimbingImageMethod = false;
    pyparams->optimizer_options.converged_force =
        params.optimizer_options.converged_force * 0.8;
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
    FILE *fileNEB = fopen(nebFilename.c_str(), "wb");
    for (long i = 0; i <= neb->numImages + 1; i++) {
      neb->path[i]->matter2con(fileNEB);
    }
    fclose(fileNEB);
    // } else
    if (status_neb == NudgedElasticBand::NEBStatus::GOOD &&
        eonc::helpers::surrogate::accuratePES(neb->path, pot)) {
      break;
    } else {
      continue;
    }
  }
  neb->printImageData();
  neb->findExtrema();
  saveData(status_neb, std::move(neb));
  return returnFiles;
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
  fileResults << magic_enum::enum_name<PotType>(
                     params.potential_options.potential)
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

  std::ofstream fileNEB(nebFilename);
  if (!fileNEB) {
    // Handle file open error
    throw std::runtime_error("Failed to open file: " + nebFilename);
  }

  for (long i = 0; i <= neb->numImages + 1; i++) {
    neb->path[i]->matter2con(nebFilename, true);
  }

  fileNEB.close();

  returnFiles.push_back("neb.dat");
  neb->printImageData(true);
}
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
        matobjs[idx].getForcesFree().array() * -1;
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
        matobjs[idx]->getForcesFree().array() * -1;
  }
  std::ostringstream oss;
  oss << targets;
  EONC_LOG_TRACE("Targets\n:{}", oss.str());
  return targets;
}
std::vector<Matter> getMidSlice(const std::vector<Matter> &matobjs) {
  // Used to get the initial data slice, endpoints and the midpoint
  std::vector<Matter> res;
  res.reserve(3);
  res.push_back(matobjs.front());
  // BUG: THIS ISN'T THE MIDDLE!!!!
  // XXX: Why does this have to be in the same order?
  // front mid back doesn't work
  // front back mid works
  res.push_back(matobjs.back());
  res.push_back(matobjs[((matobjs.size() - 2) * 2.0 / 3.0) + 1]);
  return res;
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
  Eigen::VectorXd predEnergies{Eigen::VectorXd::Zero(matobjs.size())};
  Eigen::VectorXd trueEnergies{Eigen::VectorXd::Zero(matobjs.size())};
  Eigen::VectorXd accuracy{Eigen::VectorXd::Zero(matobjs.size())};
  for (auto idx{0}; idx < predEnergies.size(); idx++) {
    predEnergies[idx] = matobjs[idx]->getPotentialEnergy();
    matobjs[idx]->setPotential(true_pot);
    trueEnergies[idx] = matobjs[idx]->getPotentialEnergy();

    accuracy[idx] = std::sqrt(predEnergies[idx] * predEnergies[idx] -
                              trueEnergies[idx] * trueEnergies[idx]);
  }
  Eigen::VectorXd difference = predEnergies - trueEnergies;
  auto mae = difference.array()
                 .abs()
                 .maxCoeff(); //.squaredNorm() / predEnergies.size();
  std::ostringstream oss;
  oss << "predicted\n"
      << predEnergies << "\ntrue\n"
      << trueEnergies << "\ndifference\n"
      << difference << "\n MAE: " << mae;
  EONC_LOG_TRACE("{}", oss.str());
  return mae < 0.05;
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
