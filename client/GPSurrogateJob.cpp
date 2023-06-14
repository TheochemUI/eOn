#include "GPSurrogateJob.h"
#include "Eigen/src/Core/Matrix.h"
#include "Job.h"
#include "Matter.h"
#include "NudgedElasticBand.h"
#include "potentials/PySurrogate/PySurrogate.h"
#include <memory>
#include <pybind11/embed.h>
#include <pybind11/eigen.h>
#include <fmt/format.h>
#include <fmt/ostream.h>

std::vector<std::string> GPSurrogateJob::run(void) {
  std::vector<std::string> returnFiles;

  // Start working
  std::string reactantFilename = helper_functions::getRelevantFile("reactant.con");
  std::string productFilename = helper_functions::getRelevantFile("product.con");

  // Clone and setup "true" params
  auto true_params = std::make_unique<Parameters>(*params);
  true_params->potential = params->true_pot;
  true_params->job = params->sub_job;

  // Get possible initial data source
  auto initial = std::make_unique<Matter>(true_params.get());
  initial->con2matter(reactantFilename);
  auto final_state = std::make_unique<Matter>(true_params.get());
  final_state->con2matter(productFilename);
  auto init_path = helper_functions::neb_paths::linearPath(*initial, *final_state, params->nebImages);
  auto init_data = helper_functions::surrogate::getMidSlice(init_path);
  // fmt::print("\nFeatures\n");
  auto features = helper_functions::surrogate::get_features(init_data);
  // fmt::print("\nAnd now the targets are\n");
  auto targets = helper_functions::surrogate::get_targets(init_data);

  // Setup a GPR Potential
  pybind11::scoped_interpreter guard{}; // Initialize the Python interpreter
  auto pot = std::make_shared<PySurrogate>(params.get());
  pot->gpmod.attr("optimize")(features, targets);
  // py::print(pot->gpmod.attr("predict")(features));
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
    // fmt::print("\nrows: {}, cols:{}\n", matobjs.size(), matobjs.front().numberOfFreeAtoms()*3);
    for (long idx{0}; idx < features.rows(); idx++){
      features.row(idx) = matobjs[idx].getPositionsFreeV();
    }
    // std::cout<<features;
    return features;
  }
  Eigen::MatrixXd get_targets(std::vector<Matter>& matobjs){
    // Always with derivatives for now
    // Energy + Derivatives for each row
    const auto nrows = matobjs.size();
    const auto ncols = (matobjs.front().numberOfFreeAtoms()*3) + 1;
    Eigen::MatrixXd targets(nrows, ncols);
    for (long idx{0}; idx < targets.rows(); idx++){
      targets.row(idx)[0] = matobjs[idx].getPotentialEnergy();
      targets.block(idx, 1, 1, ncols-1) = matobjs[idx].getForcesFreeV();
    }
    // std::cout<<targets;
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
