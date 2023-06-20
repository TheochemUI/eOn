#include "GPSurrogateJob.h"
#include "NudgedElasticBandJob.h"
#include "Potential.h"

std::vector<std::string> GPSurrogateJob::run(void) {
  std::vector<std::string> returnFiles;

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
  auto init_data = helper_functions::surrogate::getMidSlice(init_path);
  auto features = helper_functions::surrogate::get_features(init_data);
  auto targets = helper_functions::surrogate::get_targets(init_data, pot);

  // Setup a GPR Potential
  auto pypot = std::make_shared<PySurrogate>(pyparams);
  pypot->train_optimize(features, targets);
  // We only really care about the max variance, .maxCoeff()
  py::print(pypot->gpmod.attr("calculate_pred_variance")(features));
  // Make a new matter object
  // auto m1 = Matter(pypot, pyparams);
  // m1.con2matter(reactantFilename);
  // std::cout<<m1.getPotentialEnergy();
  // Start an NEB job
  // auto nebjob = NudgedElasticBandJob(pypot, pyparams);
  // nebjob.run();
  // Start an NEB run
  auto neb = std::make_unique<NudgedElasticBand>(initial, final_state, pyparams, pypot);
  auto status{neb->compute()};
  // BUG: This includes the endpoints again!!
  auto more_data = helper_functions::surrogate::get_features(neb->image);
  auto more_targets = helper_functions::surrogate::get_targets(neb->image, pot);
  auto concatFeat = helper_functions::eigen::vertCat(more_data, features);
  auto concatTargets = helper_functions::eigen::vertCat(more_targets, targets);
  while (neb->numExtrema != 2){
    auto odat = more_data;
    auto otarg = more_targets;
    more_data = helper_functions::surrogate::get_features(neb->image);
    more_targets = helper_functions::surrogate::get_targets(neb->image, pot);
    concatFeat = helper_functions::eigen::vertCat(odat, more_data);
    concatTargets = helper_functions::eigen::vertCat(otarg, more_targets);
    pypot->train_optimize(concatFeat, concatTargets);
    neb = std::make_unique<NudgedElasticBand>(initial, final_state, pyparams, pypot);
    neb->compute();
    SPDLOG_TRACE("Got num_extrema {}\n", neb->numExtrema);
  }

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
      targets.block(idx, 1, 1, ncols-1) = matobjs[idx].getForcesFree().array();
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
      targets.block(idx, 1, 1, ncols-1) = matobjs[idx]->getForcesFree().array();
    }
    SPDLOG_TRACE("Targets\n:{}", fmt::streamed(targets));
    return targets;
  }
  std::vector<Matter> getMidSlice(const std::vector<Matter>& matobjs){
    // Used to get the initial data slice, endpoints and the midpoint
    std::vector<Matter> res;
    res.reserve(3);
    res.push_back(matobjs.front());
    res.push_back(matobjs[matobjs.size()/2]);
    res.push_back(matobjs.back());
    return res;
  }
}

namespace helper_functions::eigen {
  Eigen::MatrixXd vertCat(const Eigen::MatrixXd& m1, const Eigen::MatrixXd& m2){
    assert (m1.cols() == m2.cols());
    Eigen::MatrixXd res( m1.rows() + m2.rows(), m2.cols() );
    res << m1, m2;
    return res;
  }
}
