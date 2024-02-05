#include "GPSurrogateJob.h"
#include "BaseStructures.h"
#include "Eigen/src/Core/Matrix.h"
#include "Eigen/src/Core/util/Constants.h"
#include "NudgedElasticBand.h"
#include "NudgedElasticBandJob.h"
#include "SurrogatePotential.h"
#include "helpers/Create.hpp"
#include <limits>
#include <random>

void writeDataToCSV(const std::string &filename,
                    const std::vector<double> &iterations_gp,
                    const std::vector<double> &mae_energies,
                    const std::vector<double> &true_force_norm_cis,
                    const std::vector<double> &energy_variances,
                    const std::vector<double> &rmsF_cis,
                    const std::vector<double> &maxF_cis) {
  std::ofstream csvFile(filename, std::ios::out | std::ios::trunc);
  csvFile << "Iteration,MAE_Energy,True_Force_Norm,Energy_Variance,RMSF_"
             "CI,MaxF_CI\n";
  for (size_t i = 0; i < iterations_gp.size(); ++i) {
    csvFile << fmt::format("{},{:.4e},{:.4e},{:.4e},{:.4e},{:.4e}\n",
                           iterations_gp[i], mae_energies[i],
                           true_force_norm_cis[i], energy_variances[i],
                           rmsF_cis[i], maxF_cis[i]);
  }

  csvFile.close();
}

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
  auto start = std::chrono::steady_clock::now();
  auto status_neb{neb->compute()};
  // Timer
  std::chrono::duration<double> elp_time =
      std::chrono::steady_clock::now() - start;
  SPDLOG_TRACE("Total Optimizer time for the NEB on the GP: {}s",
               elp_time.count());
  bool job_not_finished{true};
  bool train_ci{false};
  size_t n_gp{0};
  size_t n_skipat{5};
  bool retrainGPR = true;
  bool pruneOK = true;
  double pruneIncrement = 0.3;
  // Tracking variables
  std::vector<double> iterations_gp;
  std::vector<double> mae_energies;
  std::vector<double> true_force_norm_cis;
  std::vector<double> energy_variances;
  std::vector<double> rmsF_cis;
  std::vector<double> maxF_cis;
  while (job_not_finished) {
    n_gp++;
    if (n_gp > 750) {
      SPDLOG_CRITICAL("Whoops, power level of problem too high!!");
      break;
    }
    if (retrainGPR) {
      SPDLOG_TRACE("Must handle update to the GP, update number {}", n_gp);
      auto [feature, target] = helper_functions::surrogate::getNewDataPoint(
          neb->path, pot, params->gp_mindist, surpot->failedOptim);
      helper_functions::eigen::addVectorRow(features, feature);
      helper_functions::eigen::addVectorRow(targets, target);
      surpot->train_optimize(features, targets);
      for (auto &&obj : neb->path) {
        obj->setPotential(surpot);
      }
    }
    if (n_gp > 20) {
      // SPDLOG_TRACE("Using previous path due to enough data being present");
      // pyparams->gp_linear_path_always = false;
      if (n_gp % 4 == 0) {
        SPDLOG_TRACE("Halving uncertainity");
        pyparams->gp_uncertainity *= 0.5;
      }
    }
    if (!(pyparams->gp_linear_path_always) && retrainGPR) {
      SPDLOG_TRACE("Using previous path");
      neb = std::make_unique<NudgedElasticBand>(neb->path, pyparams, surpot);
    } else if (retrainGPR) {
      SPDLOG_TRACE("Using linear interpolation");
      neb = std::make_unique<NudgedElasticBand>(initial, final_state, pyparams,
                                                surpot);
    }

    start = std::chrono::steady_clock::now();
    status_neb = neb->compute();
    elp_time = std::chrono::steady_clock::now() - start;
    SPDLOG_TRACE("Total Optimizer time for the NEB on the GP: {}s",
                 elp_time.count());
    std::string nebFilename(fmt::format("neb_final_gpr_{:03d}.con", n_gp));
    returnFiles.push_back(nebFilename);
    FILE *fileNEB = fopen(nebFilename.c_str(), "wb");
    for (long i = 0; i <= neb->numImages + 1; i++) {
      neb->path[i]->matter2con(fileNEB);
    }
    fclose(fileNEB);
    if (status_neb == NudgedElasticBand::NEBStatus::GOOD) {
      if (pyparams->nebClimbingImageMethod) {
        neb->printImageDataGP(true, 0, n_gp);
        double pred_energy =
            neb->path[neb->climbingImage]->getPotentialEnergy();
        double pred_energy_variance =
            neb->path[neb->climbingImage]->getEnergyVariance();
        auto pred_forces = neb->path[neb->climbingImage]->getForcesFreeV();
        // Now get true energy, forces at CI
        neb->path[neb->climbingImage]->setPotential(pot);
        double true_energy =
            neb->path[neb->climbingImage]->getPotentialEnergy();
        // np.abs(self.energy_pred-self.energy_true)<=2.0*unc_convergence
        auto true_forces = neb->path[neb->climbingImage]->getForcesFreeV();
        double true_force_ci_norm = true_forces.norm();
        double force_ci_norm_diff = (pred_forces - true_forces).norm();
        size_t n_force_elements = true_forces.size();
        double rmsF_ci = true_force_ci_norm / std::sqrt(n_force_elements);
        // double rmsF_ci_diff = (true_forces - pred_forces).norm() / std::sqrt(n_force_elements);
        double mae_energy = abs(true_energy - pred_energy);
        double maxF_ci = abs(true_forces.maxCoeff());
        iterations_gp.push_back(n_gp);
        mae_energies.push_back(mae_energy);
        // true_force_norm_cis.push_back(true_force_ci_norm);
        true_force_norm_cis.push_back(force_ci_norm_diff);
        energy_variances.push_back(pred_energy_variance);
        rmsF_cis.push_back(rmsF_ci);
        maxF_cis.push_back(maxF_ci);
        writeDataToCSV("conv_state_gp.csv", iterations_gp, mae_energies,
                       true_force_norm_cis, energy_variances, rmsF_cis,
                       maxF_cis);

        // Display table header
        SPDLOG_TRACE("\n{:>10} {:>12} {:>18} {:>20} {:>12} {:>12} {:>12}",
                     "Iteration", "MAE Energy", "Force Diff Norm",
                     "Energy Variance", "RMSF CI", "MaxF CI", "N_GP");
        SPDLOG_TRACE(
            "---------------------------------------------------------"
            "---------------------------------------------------------");
        // Display each row
        for (size_t idx = 0; idx < mae_energies.size(); ++idx) {
          SPDLOG_TRACE(
              "{:>10} {:>12.4e} {:>18.4e} {:>20.4e} {:>12.4e} {:>12.4e} {:>10}",
              idx + 1, mae_energies[idx], true_force_norm_cis[idx],
              energy_variances[idx], rmsF_cis[idx], maxF_cis[idx],
              iterations_gp[idx]);
        }

        // 0.0003 Eh/Bohr is around 0.01543 eV/A
        // 0.0005 Eh/Bohr is around 0.02571 eV/A
        if ((mae_energies.back() < 0.01543) || (true_force_norm_cis.back() < 0.02571)) {
          SPDLOG_INFO("Converged due to low force and energy differences on "
                      "true surface at the CI");
          break;
        }
      } else {
        if (iterations_gp.size() > 2) {
          if ((n_gp > iterations_gp.front() + 3) && pruneOK) {
            SPDLOG_INFO("Pruning since {} and we first converged at {}", n_gp,
                        iterations_gp.front());
            helper_functions::surrogate::pruneHighForceData(
                features, targets, n_gp - iterations_gp.front());
            pruneOK = false;
          }
        }
        pyparams->nebClimbingImageMethod = true;
        pyparams->nebClimbingImageConvergedOnly = true;
        pyparams->gp_linear_path_always = false;
        retrainGPR = true;
      }
    }
  }
  // Just to get better images, will use the true potential
  pot->m_log->info("Switching to true potential for final extreum");
  for (auto &&obj : neb->path) {
    obj->setPotential(pot);
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
  // SPDLOG_TRACE("Features\n:{}", fmt::streamed(features));
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
    targets(idx, Eigen::seqN(1, ncols - 1)) =
        matobjs[idx].getForcesFreeV() * -1; // gradients
  }
  // SPDLOG_TRACE("Targets\n:{}", fmt::streamed(targets));
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
void addCI(Eigen::MatrixXd &features, Eigen::MatrixXd &targets,
           NudgedElasticBand *neb, double true_energy,
           Eigen::VectorXd true_forces) {
  features.conservativeResize(features.rows() + 1, Eigen::NoChange);
  features.row(features.rows() - 1) =
      neb->path[neb->climbingImage]->getPositionsV();
  targets.conservativeResize(targets.rows() + 1, Eigen::NoChange);
  targets(targets.rows() - 1, 0) = true_energy;
  targets.block(targets.rows() - 1, 1, 1, true_forces.size()) =
      -1 * true_forces;
}
std::pair<Eigen::VectorXd, Eigen::VectorXd>
getNewDataPoint(const std::vector<std::shared_ptr<Matter>> &matobjs,
                std::shared_ptr<Potential> true_pot,
                double min_distance_threshold, bool optfail) {

  auto [maxUnc, maxIndex] = getMaxUncertainty(matobjs);
  Matter candidate{*matobjs[maxIndex + 1]};

  // Check if the point is too close to existing points
  bool tooClose = false;
  double min_cdistext = std::numeric_limits<double>::infinity();
  for (const auto &existing_obj : matobjs) {
    double cdistext = candidate.distanceTo(*existing_obj);
    if (cdistext < 1e-8) { // the point itself
      continue;
    }
    min_cdistext = std::min(min_cdistext, cdistext);
  }
  SPDLOG_INFO("Got minimal distance {}, compared to {}", min_cdistext,
              min_distance_threshold);
  if (min_cdistext > min_distance_threshold) {
    tooClose = true;
  }

  // If too close or optimizer failed, select a random point along the path
  if (tooClose || optfail) {
    std::random_device rd;
    std::mt19937 gen(rd());
    // No endpoints
    std::uniform_int_distribution<> distrib(1, matobjs.size() - 2);
    int randomIndex = distrib(gen);
    candidate = *matobjs[randomIndex];
    SPDLOG_INFO("Optimizer failed or points too close, picking random index {}",
                randomIndex);
  }
  // SPDLOG_INFO("Got {}", fmt::streamed(candidate.getPositionsFreeV()));

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

bool pruneHighForceData(Eigen::MatrixXd &features, Eigen::MatrixXd &targets,
                        int fixedRowsToKeep) {
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
