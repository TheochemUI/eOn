#include "GPSurrogateJob.h"
#include "BaseStructures.h"
#include "NudgedElasticBand.h"
#include "NudgedElasticBandJob.h"
#include "SurrogatePotential.h"
#include "helpers/Create.hpp"

std::vector<std::string> GPSurrogateJob::run(void) {
  // Start working
  std::string reactantFilename =
      helper_functions::getRelevantFile("reactant.con");
  std::string productFilename =
      helper_functions::getRelevantFile("product.con");

  // Clone and setup "true" params
  auto true_params = std::make_shared<Parameters>(*params);
  true_params->job = params->sub_job;
  auto true_job =
      helper_functions::makeJob(std::make_unique<Parameters>(*true_params));
  auto pyparams = std::make_shared<Parameters>(*params);
  pyparams->potential = params->surrogatePotential;
  // pyparams->nebImages = params->nebImages + 2;

  // Get possible initial data source
  auto initial = std::make_shared<Matter>(pot, true_params);
  initial->con2matter(reactantFilename);
  auto final_state = std::make_shared<Matter>(pot, true_params);
  final_state->con2matter(productFilename);
  auto init_path = helper_functions::neb_paths::linearPath(
      *initial, *final_state, params->nebImages);
  auto init_data = helper_functions::surrogate::getMidSlice(init_path);
  auto features = helper_functions::surrogate::get_features(init_data);
  SPDLOG_TRACE("Potential is {}",
               helper_functions::getPotentialName(pot->getType()));
  auto targets = helper_functions::surrogate::get_targets(init_data, pot);

  // Setup a GPR Potential
  auto surpot = helpers::create::makeSurrogatePotential(
      params->surrogatePotential, params, initial);
  surpot->train_optimize(features, targets);
  pyparams->nebClimbingImageMethod = false;
  pyparams->nebClimbingImageConvergedOnly = false;
  auto neb = std::make_unique<NudgedElasticBand>(initial, final_state, pyparams,
                                                 surpot);
  auto status_neb{neb->compute()};
  // helper_functions::surrogate::accuratePES(neb->path, pot, params->gp_accuracy);
  bool job_not_finished{true};
  size_t n_gp{0};
  int nrow{5};
  pyparams->gp_uncertainity = 0.1; // Start out high
  bool retrainGPR = true;
  while (job_not_finished) { // outer loop?
    n_gp++;
    if (n_gp > 750) {
      SPDLOG_CRITICAL("Whoops, power level of problem too high!!");
      break;
    }
    // if (status_neb == NudgedElasticBand::NEBStatus::MAX_UNCERTAINITY ||
    // status_neb == NudgedElasticBand::NEBStatus::BAD_MAX_ITERATIONS) {
    if (retrainGPR) {
        SPDLOG_TRACE("Must handle update to the GP, update number {}", n_gp);
        auto [feature, target] = helper_functions::surrogate::getNewDataPoint(neb->path, pot);
        helper_functions::eigen::addVectorRow(features, feature);
        helper_functions::eigen::addVectorRow(targets, target);
        if (n_gp % 3 == 2) {
          SPDLOG_TRACE("Trying to prune");
          if (helper_functions::surrogate::pruneHighForceData(features, targets,
                                                              nrow)) {
            pyparams->gp_uncertainity /= 2;
            nrow *= 2;
            nrow = min(nrow, 36000);
            pyparams->gp_uncertainity = max(pyparams->gp_uncertainity, 0.05);
          }
        }
        surpot->train_optimize(features, targets);
        for (auto &&obj : neb->path) {
            obj->setPotential(surpot);
        }
    }
    if (!(pyparams->gp_linear_path_always) && retrainGPR) {
        SPDLOG_TRACE("Using previous path");
        neb = std::make_unique<NudgedElasticBand>(neb->path, pyparams, surpot);
    } else if (retrainGPR) {
        SPDLOG_TRACE("Using linear interpolation");
        neb = std::make_unique<NudgedElasticBand>(initial, final_state, pyparams, surpot);
    }
    status_neb = neb->compute();
    // helper_functions::surrogate::accuratePES(neb->path, pot, params->gp_accuracy);

    std::string nebFilename(fmt::format("neb_final_gpr_{:03d}.con", n_gp));
    returnFiles.push_back(nebFilename);
    FILE *fileNEB = fopen(nebFilename.c_str(), "wb");
    for (long i = 0; i <= neb->numImages + 1; i++) {
      neb->path[i]->matter2con(fileNEB);
    }
    fclose(fileNEB);
    // } else

    if (status_neb == NudgedElasticBand::NEBStatus::GOOD &&
        helper_functions::surrogate::accuratePES(neb->path, pot, params->gp_accuracy)) {
        if (!pyparams->nebClimbingImageMethod) {
            SPDLOG_TRACE("Turning on CI");
            pyparams->nebClimbingImageMethod = true;
            pyparams->nebClimbingImageConvergedOnly = true;
            pyparams->gp_uncertainity = 0.01;
            pyparams->gp_linear_path_always = false;
            retrainGPR = false; // Do not retrain after first convergence
            pyparams->nebImages = params->nebImages;
            neb = std::make_unique<NudgedElasticBand>(neb->path, pyparams, surpot);
            continue;
        } else {
            break; // Exit the loop if converged with climbing image method active
        }
    } else {
        retrainGPR = true; // Set flag to retrain GPR if not converged
    }
  }
  neb->printImageData();
  neb->findExtrema();
  saveData(status_neb, std::move(neb));
  return returnFiles;
}

void GPSurrogateJob::saveData(NudgedElasticBand::NEBStatus status,
                              std::unique_ptr<NudgedElasticBand> neb) {
  std::string resultsFilename = "results_gpr.dat";
  returnFiles.push_back(resultsFilename);

  std::ofstream fileResults(resultsFilename);
  if (!fileResults) {
    // Handle file open error
    throw std::runtime_error("Failed to open file: " + resultsFilename);
  }

  fileResults << static_cast<int>(status) << " termination_reason\n";
  fileResults << helper_functions::getPotentialName(params->potential)
              << " potential_type\n";
  fileResults << fmt::format("{:.6f} energy_reference\n",
                             neb->path[0]->getPotentialEnergy());
  fileResults << neb->numImages << " number_of_images\n";

  for (long i = 0; i <= neb->numImages + 1; i++) {
    fileResults << fmt::format("{:.6f} image{}_energy\n",
                               neb->path[i]->getPotentialEnergy() -
                                   neb->path[0]->getPotentialEnergy(),
                               i);
    fileResults << fmt::format("{:.6f} image{}_force\n",
                               neb->path[i]->getForces().norm(), i);
    fileResults << fmt::format("{:.6f} image{}_projected_force\n",
                               neb->projectedForce[i]->norm(), i);
  }

  fileResults << neb->numExtrema << " number_of_extrema\n";
  for (long i = 0; i < neb->numExtrema; i++) {
    fileResults << fmt::format("{:.6f} extremum{}_position\n",
                               neb->extremumPosition[i], i);
    fileResults << fmt::format("{:.6f} extremum{}_energy\n",
                               neb->extremumEnergy[i], i);
  }

  fileResults.close();

  std::string nebFilename = "neb_gpr.con";
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

  returnFiles.push_back("neb_gpr.dat");
  neb->printImageData(true);
}
namespace helper_functions::surrogate {
Eigen::MatrixXd get_features(const std::vector<Matter> &matobjs) {
  // Calculate dimensions
  Eigen::MatrixXd features(matobjs.size(),
                           matobjs.front().numberOfFreeAtoms() * 3);
  SPDLOG_TRACE("rows: {}, cols:{}", matobjs.size(),
               matobjs.front().numberOfFreeAtoms() * 3);
  for (long idx{0}; idx < features.rows(); idx++) {
    features.row(idx) = matobjs[idx].getPositionsFreeV();
  }
  SPDLOG_TRACE("Features\n:{}", fmt::streamed(features));
  return features;
}
Eigen::MatrixXd
get_features(const std::vector<std::shared_ptr<Matter>> &matobjs) {
  // Calculate dimensions
  Eigen::MatrixXd features(matobjs.size(),
                           matobjs.front()->numberOfFreeAtoms() * 3);
  SPDLOG_TRACE("rows: {}, cols:{}\n", matobjs.size(),
               matobjs.front()->numberOfFreeAtoms() * 3);
  for (long idx{0}; idx < features.rows(); idx++) {
    features.row(idx) = matobjs[idx]->getPositionsFreeV();
  }
  SPDLOG_TRACE("Features\n:{}", fmt::streamed(features));
  return features;
}
Eigen::MatrixXd get_targets(std::vector<Matter> &matobjs,
                            std::shared_ptr<Potential> true_pot) {
  // Always with derivatives for now
  // Energy + gradients for each row
  const auto nrows = matobjs.size();
  const auto ncols = (matobjs.front().numberOfFreeAtoms() * 3) + 1;
  Eigen::MatrixXd targets(nrows, ncols);
  for (long idx{0}; idx < targets.rows(); idx++) {
    matobjs[idx].setPotential(true_pot);
    targets.row(idx)[0] = matobjs[idx].getPotentialEnergy();
    targets(idx, Eigen::seqN(1, ncols-1)) =
        matobjs[idx].getForcesFreeV() * -1; // gradients
  }
  SPDLOG_TRACE("Targets\n:{}", fmt::streamed(targets));
  return targets;
}
Eigen::MatrixXd get_targets(std::vector<std::shared_ptr<Matter>> &matobjs,
                            std::shared_ptr<Potential> true_pot) {
  const auto nrows = matobjs.size();
  const auto ncols = (matobjs.front()->numberOfFreeAtoms() * 3) + 1;
  Eigen::MatrixXd targets(nrows, ncols);
  for (long idx{0}; idx < targets.rows(); idx++) {
    matobjs[idx]->setPotential(true_pot);
    targets.row(idx)[0] = matobjs[idx]->getPotentialEnergy();
    targets.block(idx, 1, 1, ncols - 1) =
        matobjs[idx]->getForcesFree().array() * -1;
  }
  SPDLOG_TRACE("Targets\n:{}", fmt::streamed(targets));
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
  // SPDLOG_TRACE("Generated Target:\n{}", fmt::streamed(target));
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
  // SPDLOG_TRACE("Uncertainity along path is {}\nmax_index: {}, maxVal: {}",
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
                 std::shared_ptr<Potential> true_pot, double max_accuracy) {
  Eigen::VectorXd predEnergies(matobjs.size());
  Eigen::VectorXd trueEnergies(matobjs.size());

  for (int idx = 0; idx < matobjs.size(); idx++) {
    predEnergies[idx] = matobjs[idx]->getPotentialEnergy();
    matobjs[idx]->setPotential(true_pot);
    trueEnergies[idx] = matobjs[idx]->getPotentialEnergy();
  }

  Eigen::VectorXd difference = predEnergies - trueEnergies;
  double mse = difference.squaredNorm() / matobjs.size();
  SPDLOG_TRACE("predicted\n{}\ntrue\n{}\ndifference\n{}\n RMSE: {}\n",
               fmt::streamed(predEnergies), fmt::streamed(trueEnergies),
               fmt::streamed(difference), sqrt(mse));
  return sqrt(mse) < max_accuracy;
}

bool pruneHighForceData(Eigen::MatrixXd &features, Eigen::MatrixXd &targets, int fixedRowsToKeep) {
  SPDLOG_TRACE("To keep: {}", fixedRowsToKeep);
    assert(features.rows() == targets.rows());

    // If the number of rows is less than or equal to the fixed size, return early
    if (features.rows() <= fixedRowsToKeep) {
        return false;
    }

    std::vector<std::pair<double, int>> forceNorms;
    for (int i = 0; i < targets.rows(); ++i) {
        double forceNorm = targets.row(i).tail(targets.cols() - 1).norm();
        forceNorms.emplace_back(forceNorm, i);
    }

    // TODO: Use the perpendicular forces
    // Sort the pairs by force norm in ascending order
    std::sort(forceNorms.begin(), forceNorms.end());

    // Keep only the fixed number of rows with the lowest force norms
    std::vector<int> rowsToKeep;
    for (int i = 0; i < fixedRowsToKeep; ++i) {
        rowsToKeep.push_back(forceNorms[i].second);
        SPDLOG_TRACE("Force norms {} kept", forceNorms[i].second);
    }

    // Create new matrices for pruned features and targets
    Eigen::MatrixXd newFeatures(fixedRowsToKeep, features.cols());
    Eigen::MatrixXd newTargets(fixedRowsToKeep, targets.cols());

    for (size_t i = 0; i < rowsToKeep.size(); ++i) {
        newFeatures.row(i) = features.row(rowsToKeep[i]);
        newTargets.row(i) = targets.row(rowsToKeep[i]);
    }
    fmt::print("Post pruning {}\n", rowsToKeep.size());

    // Replace the old matrices with the new pruned ones
    features = std::move(newFeatures);
    targets = std::move(newTargets);
    return true;
}

} // namespace helper_functions::surrogate

namespace helper_functions::eigen {
Eigen::MatrixXd vertCat(const Eigen::MatrixXd &m1, const Eigen::MatrixXd &m2) {
  assert(m1.cols() == m2.cols());
  Eigen::MatrixXd res(m1.rows() + m2.rows(), m2.cols());
  res << m1, m2;
  return res;
}
void addVectorRow(Eigen::MatrixXd &data, const Eigen::VectorXd &newrow) {
  assert(data.cols() == newrow.size());
  data.conservativeResize(data.rows() + 1, data.cols());
  data.row(data.rows() - 1) = newrow;
}
} // namespace helper_functions::eigen
