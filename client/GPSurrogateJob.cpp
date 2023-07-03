#include "GPSurrogateJob.h"
#include "NudgedElasticBand.h"
#include "NudgedElasticBandJob.h"
#include "Potential.h"

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
  pyparams->potential = PotType::PYSURROGATE;

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
  auto pypot = std::make_shared<PySurrogate>(pyparams);
  pypot->train_optimize(features, targets);
  auto neb = std::make_unique<NudgedElasticBand>(initial, final_state, pyparams,
                                                 pypot);
  auto status_neb{neb->compute()};
  bool job_not_finished{true};
  size_t n_gp{0};
  double unc_conv{pyparams->gp_uncertainity};
  while (job_not_finished) { // outer loop?
    n_gp++;
    if (n_gp > 750) {
      SPDLOG_CRITICAL("Whoops, power level of problem too high!!");
      break;
    }
    // if (status_neb == NudgedElasticBand::NEBStatus::MAX_UNCERTAINITY ||
    // status_neb == NudgedElasticBand::NEBStatus::BAD_MAX_ITERATIONS) {
    SPDLOG_TRACE("Must handle update to the GP, update number {}", n_gp);
    auto [maxUnc, maxIndex] =
        helper_functions::surrogate::getMaxUncertainty(neb->path);
    // pyparams->gp_uncertainity = ( maxUnc + unc_conv ) / n_gp++;
    SPDLOG_TRACE("New allowed uncertainity is {}", pyparams->gp_uncertainity);
    auto [feature, target] =
        helper_functions::surrogate::getNewDataPoint(neb->path, pot);
    helper_functions::eigen::addVectorRow(features, feature);
    helper_functions::eigen::addVectorRow(targets, target);
    pypot->train_optimize(features, targets);
    pyparams->nebClimbingImageMethod = false;
    pyparams->optConvergedForce = params->optConvergedForce * 0.8;
    for (auto &&obj : neb->path) {
      obj->setPotential(pypot);
    }
    neb = std::make_unique<NudgedElasticBand>(neb->path, pyparams, pypot);
    status_neb = neb->compute();

    std::string nebFilename(fmt::format("neb_final_gpr_{:03d}.con", n_gp));
    returnFiles.push_back(nebFilename);
    FILE *fileNEB = fopen(nebFilename.c_str(), "wb");
    for (long i = 0; i <= neb->numImages + 1; i++) {
      neb->path[i]->matter2con(fileNEB);
    }
    fclose(fileNEB);
    // } else
    if (status_neb == NudgedElasticBand::NEBStatus::GOOD) {
      SPDLOG_TRACE("Regular NEB converged on GP surface");
      auto pyparam_ci = std::make_shared<Parameters>(*pyparams);
      pyparam_ci->nebClimbingImageMethod = true;
      pyparam_ci->writeMovies = true;
      SPDLOG_TRACE("Starting a CI NEB run");

      for (auto &&obj : neb->path) {
        obj->setPotential(pypot);
      }
      neb = std::make_unique<NudgedElasticBand>(neb->path, pyparam_ci, pypot);
      status_neb = neb->compute();
      if (status_neb == NudgedElasticBand::NEBStatus::GOOD) {
        job_not_finished =
            !helper_functions::surrogate::accuratePES(neb->path, pot);
      } else {
        pyparams->optMethod = "lbfgs";
      }
    }
  }
  neb->printImageData();
  neb->findExtrema();
  saveData(status_neb, pot, std::move(neb));
  return returnFiles;
}

void GPSurrogateJob::saveData(NudgedElasticBand::NEBStatus status,
                              std::shared_ptr<Potential> true_pot,
                              std::unique_ptr<NudgedElasticBand> neb) {
  FILE *fileResults, *fileNEB;

  std::string resultsFilename("results.dat");
  returnFiles.push_back(resultsFilename);
  fileResults = fopen(resultsFilename.c_str(), "wb");
  // for (auto&& nebo : neb->path){
  //     nebo->setPotential(true_pot);
  // }
  fprintf(fileResults, "%d termination_reason\n", static_cast<int>(status));
  fprintf(fileResults, "%s potential_type\n",
          helper_functions::getPotentialName(params->potential).c_str());
  // fprintf(fileResults, "%ld total_force_calls\n", Potential::fcalls);
  // fprintf(fileResults, "%ld force_calls_neb\n", fCallsNEB);
  fprintf(fileResults, "%f energy_reference\n",
          neb->path[0]->getPotentialEnergy());
  fprintf(fileResults, "%li number_of_images\n", neb->numImages);
  for (long i = 0; i <= neb->numImages + 1; i++) {
    fprintf(fileResults, "%f image%li_energy\n",
            neb->path[i]->getPotentialEnergy() -
                neb->path[0]->getPotentialEnergy(),
            i);
    fprintf(fileResults, "%f image%li_force\n",
            neb->path[i]->getForces().norm(), i);
    fprintf(fileResults, "%f image%li_projected_force\n",
            neb->projectedForce[i]->norm(), i);
  }
  fprintf(fileResults, "%li number_of_extrema\n", neb->numExtrema);
  for (long i = 0; i < neb->numExtrema; i++) {
    fprintf(fileResults, "%f extremum%li_position\n", neb->extremumPosition[i],
            i);
    fprintf(fileResults, "%f extremum%li_energy\n", neb->extremumEnergy[i], i);
  }

  fclose(fileResults);

  std::string nebFilename(fmt::format("neb.con"));
  returnFiles.push_back(nebFilename);
  fileNEB = fopen(nebFilename.c_str(), "wb");
  for (long i = 0; i <= neb->numImages + 1; i++) {
    neb->path[i]->matter2con(fileNEB);
  }
  fclose(fileNEB);

  returnFiles.push_back("neb.dat");
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
  // Energy + Derivatives for each row
  const auto nrows = matobjs.size();
  const auto ncols = (matobjs.front().numberOfFreeAtoms() * 3) + 1;
  Eigen::MatrixXd targets(nrows, ncols);
  for (long idx{0}; idx < targets.rows(); idx++) {
    matobjs[idx].setPotential(true_pot);
    targets.row(idx)[0] = matobjs[idx].getPotentialEnergy();
    targets.block(idx, 1, 1, ncols - 1) =
        matobjs[idx].getForcesFree().array() * -1;
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
  return std::make_pair(maxUnc, maxIndex);
}
std::pair<Eigen::VectorXd, Eigen::VectorXd>
getNewDataPoint(const std::vector<std::shared_ptr<Matter>> &matobjs,
                std::shared_ptr<Potential> true_pot) {
  auto [maxUnc, maxIndex] = getMaxUncertainty(matobjs);
  return std::make_pair<Eigen::VectorXd, Eigen::VectorXd>(
      matobjs[maxIndex + 1]->getPositionsFreeV(),
      make_target(*matobjs[maxIndex], true_pot));
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
  double largest_error = accuracy.maxCoeff();
  SPDLOG_TRACE(
      "predicted\n{}true\n{}accuracy\n{}\nlargest_error {} will return {}",
      fmt::streamed(predEnergies), fmt::streamed(trueEnergies),
      fmt::streamed(accuracy), largest_error, largest_error < 0.05);
  return largest_error < 0.05;
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
