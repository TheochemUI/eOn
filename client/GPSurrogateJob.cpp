#include "GPSurrogateJob.h"
#include "NudgedElasticBand.h"
#include "NudgedElasticBandJob.h"
#include "Potential.h"

std::vector<std::string> GPSurrogateJob::run(void) {
  // Start working
  std::string reactantFilename = helper_functions::getRelevantFile("reactant.con");
  std::string productFilename = helper_functions::getRelevantFile("product.con");

  // Clone and setup "true" params
  auto true_params = std::make_shared<Parameters>(*params);
  true_params->job = params->sub_job;
  auto true_job = helper_functions::makeJob(std::make_unique<Parameters>(*true_params));
  auto pyparams = std::make_shared<Parameters>(*params);
  pyparams->potential = PotType::PYSURROGATE;

  // Get possible initial data source
  auto initial = std::make_shared<Matter>(pot, true_params);
  initial->con2matter(reactantFilename);
  auto final_state = std::make_shared<Matter>(pot, true_params);
  final_state->con2matter(productFilename);
  auto init_path = helper_functions::neb_paths::linearPath(*initial, *final_state, params->nebImages);
  for (std::size_t i = 0; i < init_path.size(); ++i) {
    auto &&mobj = init_path[i];
    // SPDLOG_TRACE("Index: {}, Positions Free:\n {}", i,
    //              fmt::streamed(mobj.getPositionsFree()));
  }
  auto init_data = helper_functions::surrogate::getMidSlice(init_path);
  auto features = helper_functions::surrogate::get_features(init_data);
  SPDLOG_TRACE("Potential is {}", helper_functions::getPotentialName(pot->getType()));
  auto targets = helper_functions::surrogate::get_targets(init_data, pot);

  // Setup a GPR Potential
  auto pypot = std::make_shared<PySurrogate>(pyparams);
  pypot->train_optimize(features, targets);
  // We only really care about the max variance, .maxCoeff()
  // py::print(pypot->gpmod.attr("calculate_pred_variance")(features));
  // Make a new matter object
  // auto m1 = Matter(pypot, pyparams);
  // m1.con2matter(reactantFilename);
  // std::cout<<m1.getPotentialEnergy();
  // Start an NEB job
  // auto nebjob = NudgedElasticBandJob(pypot, pyparams);
  // nebjob.run();
  // Start an NEB run
  auto neb = std::make_unique<NudgedElasticBand>(initial, final_state, pyparams, pypot);
  auto status_neb{neb->compute()};
  bool job_not_finished{true};
  size_t n_gp{0};
  while (job_not_finished) {
    n_gp++;
    if (status_neb == NudgedElasticBand::NEBStatus::MAX_UNCERTAINITY) {
      SPDLOG_TRACE("Must handle update to the GP, update number {}", n_gp);
      auto [feature, target] =
          helper_functions::surrogate::getNewDataPoint(neb->image);
      helper_functions::eigen::addVectorRow(features, feature);
      // SPDLOG_TRACE("New Features:\n {}", fmt::streamed(features));
      helper_functions::eigen::addVectorRow(targets, target);
      // SPDLOG_TRACE("New Targets:\n {}", fmt::streamed(targets));
      pypot->train_optimize(features, targets);
      // if ( status==NudgedElasticBand::NEBStatus::BAD_MAX_ITERATIONS && n_gp >
      // 5){
      //   // std::exit(1);
      //   // pyparams->optMaxMove += 0.05;
      //   // pyparams->optMaxIterations += 50;
      //   // status = neb->compute();
      // } else {
      //   neb = std::make_unique<NudgedElasticBand>(initial, final_state,
      //   pyparams, pypot); status = neb->compute();
      // }
      neb = std::make_unique<NudgedElasticBand>(initial, final_state, pyparams,
                                                pypot);
      status_neb = neb->compute();
    }
    if (status_neb == NudgedElasticBand::NEBStatus::GOOD) {
      job_not_finished = false;
    }
  }
  saveData(status_neb, pot, std::move(neb));
  // BUG: This includes the endpoints again!!
  // auto more_data = helper_functions::surrogate::get_features(neb->image);
  // auto more_targets = helper_functions::surrogate::get_targets(neb->image, pot);
  // auto concatFeat = helper_functions::eigen::vertCat(more_data, features);
  // auto concatTargets = helper_functions::eigen::vertCat(more_targets, targets);
  // while (neb->numExtrema != 2){
  //   auto odat = more_data;
  //   auto otarg = more_targets;
  //   more_data = helper_functions::surrogate::get_features(neb->image);
  //   more_targets = helper_functions::surrogate::get_targets(neb->image, pot);
  //   concatFeat = helper_functions::eigen::vertCat(odat, more_data);
  //   concatTargets = helper_functions::eigen::vertCat(otarg, more_targets);
  //   pypot->train_optimize(concatFeat, concatTargets);
  //   neb = std::make_unique<NudgedElasticBand>(initial, final_state, pyparams, pypot);
  //   neb->compute();
  //   SPDLOG_TRACE("Got num_extrema {}\n", neb->numExtrema);
  // }

  // py::print(pypot->gpmod.attr("predict")(features));
  // initial->setPotential(pypot);
  // std::cout<<(initial->getForces());
  // auto [energy, forces] = pot->get_ef(initial->getPositions(), initial->getAtomicNrs(), initial->getCell());
  // auto execString =
  //   fmt::format("Energy:\n{energy:}\nForces:\n{forces:}",
  //                 fmt::arg("energy", energy),
  //                 fmt::arg("forces", fmt::streamed(forces)));
  // std::cout << execString << "\n";

  // auto sub_job = helper_functions::makeJob(std::move(true_params));
  // sub_job->run();
  return returnFiles;
}

void GPSurrogateJob::saveData(NudgedElasticBand::NEBStatus status,
                              std::shared_ptr<Potential> true_pot,
                              std::unique_ptr<NudgedElasticBand> neb) {
  FILE *fileResults, *fileNEB;

  std::string resultsFilename("results.dat");
  returnFiles.push_back(resultsFilename);
  fileResults = fopen(resultsFilename.c_str(), "wb");
  for (auto&& nebo : neb->image){
      nebo->setPotential(true_pot);
  }
  fprintf(fileResults, "%d termination_reason\n", static_cast<int>(status));
  fprintf(fileResults, "%s potential_type\n",
          helper_functions::getPotentialName(params->potential).c_str());
  // fprintf(fileResults, "%ld total_force_calls\n", Potential::fcalls);
  // fprintf(fileResults, "%ld force_calls_neb\n", fCallsNEB);
  fprintf(fileResults, "%f energy_reference\n",
          neb->image[0]->getPotentialEnergy());
  fprintf(fileResults, "%li number_of_images\n", neb->images);
  for (long i = 0; i <= neb->images + 1; i++) {
    fprintf(fileResults, "%f image%li_energy\n",
            neb->image[i]->getPotentialEnergy() -
                neb->image[0]->getPotentialEnergy(),
            i);
    fprintf(fileResults, "%f image%li_force\n",
            neb->image[i]->getForces().norm(), i);
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

  std::string nebFilename("neb.con");
  returnFiles.push_back(nebFilename);
  fileNEB = fopen(nebFilename.c_str(), "wb");
  for (long i = 0; i <= neb->images + 1; i++) {
    neb->image[i]->matter2con(fileNEB);
  }
  fclose(fileNEB);

  returnFiles.push_back("neb.dat");
  neb->printImageData(true);
}

namespace helper_functions::surrogate {
  Eigen::MatrixXd get_features(const std::vector<Matter>& matobjs){
    // Calculate dimensions
    Eigen::MatrixXd features(matobjs.size(), matobjs.front().numberOfFreeAtoms()*3);
    SPDLOG_TRACE("rows: {}, cols:{}", matobjs.size(), matobjs.front().numberOfFreeAtoms()*3);
    for (long idx{0}; idx < features.rows(); idx++){
      features.row(idx) = matobjs[idx].getPositionsFreeV();
    }
    SPDLOG_TRACE("Features\n:{}", fmt::streamed(features));
    return features;
  }
  Eigen::MatrixXd get_features(const std::vector<std::shared_ptr<Matter>>& matobjs){
    // Calculate dimensions
    Eigen::MatrixXd features(matobjs.size(), matobjs.front()->numberOfFreeAtoms()*3);
    SPDLOG_TRACE("rows: {}, cols:{}\n", matobjs.size(), matobjs.front()->numberOfFreeAtoms()*3);
    for (long idx{0}; idx < features.rows(); idx++){
      features.row(idx) = matobjs[idx]->getPositionsFreeV();
    }
    SPDLOG_TRACE("Features\n:{}", fmt::streamed(features));
    return features;
  }
  Eigen::MatrixXd get_targets(std::vector<Matter>& matobjs, std::shared_ptr<Potential> true_pot){
    // Always with derivatives for now
    // Energy + Derivatives for each row
    const auto nrows = matobjs.size();
    const auto ncols = (matobjs.front().numberOfFreeAtoms()*3) + 1;
    Eigen::MatrixXd targets(nrows, ncols);
    for (long idx{0}; idx < targets.rows(); idx++){
      matobjs[idx].setPotential(true_pot);
      targets.row(idx)[0] = matobjs[idx].getPotentialEnergy();
      targets.block(idx, 1, 1, ncols-1) = matobjs[idx].getForcesFree().array() * -1;
    }
    SPDLOG_TRACE("Targets\n:{}", fmt::streamed(targets));
    return targets;
  }
  Eigen::MatrixXd get_targets(std::vector<std::shared_ptr<Matter>>& matobjs, std::shared_ptr<Potential> true_pot){
    const auto nrows = matobjs.size();
    const auto ncols = (matobjs.front()->numberOfFreeAtoms()*3) + 1;
    Eigen::MatrixXd targets(nrows, ncols);
    for (long idx{0}; idx < targets.rows(); idx++){
      matobjs[idx]->setPotential(true_pot);
      targets.row(idx)[0] = matobjs[idx]->getPotentialEnergy();
      targets.block(idx, 1, 1, ncols-1) = matobjs[idx]->getForcesFree().array() * -1;
    }
    SPDLOG_TRACE("Targets\n:{}", fmt::streamed(targets));
    return targets;
  }
  std::vector<Matter> getMidSlice(const std::vector<Matter>& matobjs){
    // Used to get the initial data slice, endpoints and the midpoint
    std::vector<Matter> res;
    res.reserve(3);
    res.push_back(matobjs.front());
    // BUG: THIS ISN'T THE MIDDLE!!!!
    // XXX: Why does this have to be in the same order?
    // front mid back doesn't work
    // front back mid works
    res.push_back(matobjs.back());
    res.push_back(matobjs[(( matobjs.size() - 2 )*2.0/3.0)+1]);
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
  std::pair<Eigen::VectorXd, Eigen::VectorXd>
  getNewDataPoint(const std::vector<std::shared_ptr<Matter>> &matobjs,
                  std::shared_ptr<Potential> true_pot) {
    // TODO: Refactor
    // This function takes the path, finds the "best" new data point to be
    // evaluated and returns the features and targets for that particular point
    // NOTE: This assumes that the Surrogate potential is in use
    Eigen::VectorXd pathUncertainity{Eigen::VectorXd::Zero(matobjs.size())};
    for (auto idx {0}; idx < matobjs.size(); idx++) {
      pathUncertainity[idx] = matobjs[idx]->getEnergyVariance();
    }
    int maxIndex = 0;
    for (int idx = 0; idx < pathUncertainity.size(); idx++) {
      if (pathUncertainity[idx] == pathUncertainity.maxCoeff()) {
        maxIndex = idx;
        break;
      }
    }
    SPDLOG_TRACE("Max uncertainity is {}\n from {} at {}", pathUncertainity.maxCoeff(),
                 fmt::streamed(pathUncertainity), maxIndex);
    // SPDLOG_TRACE("Generated Feature:\n{}", fmt::streamed(matobjs[maxIndex]->getPositionsFreeV()));
    return std::make_pair<Eigen::VectorXd, Eigen::VectorXd>(
        matobjs[maxIndex]->getPositionsFreeV(),
        make_target(*matobjs[maxIndex], true_pot));
  }
}

namespace helper_functions::eigen {
  Eigen::MatrixXd vertCat(const Eigen::MatrixXd& m1, const Eigen::MatrixXd& m2){
    assert (m1.cols() == m2.cols());
    Eigen::MatrixXd res( m1.rows() + m2.rows(), m2.cols() );
    res << m1, m2;
    return res;
  }
  void addVectorRow(Eigen::MatrixXd& data, const Eigen::VectorXd& newrow){
   assert (data.cols() == newrow.size());
   data.conservativeResize(data.rows()+1, data.cols());
   data.row(data.rows() - 1) = newrow;
  }
}
