#include "BaseStructures.h"
#include "GPSurrogateJob.h"
#include "Matter.h"
#include "NudgedElasticBand.h"
#include "Parameters.h"
#include "helpers/Create.hpp"
#include "potentials/GPRPotential/GPRPotential.h"
#include <cstdlib>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <memory>
#include <pybind11/embed.h>

using namespace std::string_literals; // For ""s

int main(void) {
  auto params = std::make_shared<Parameters>();
  auto true_pot = helper_functions::makePotential(PotType::CUH2, params);
  // Make features
  auto initial = std::make_shared<Matter>(true_pot, params);
  initial->con2matter("reactant.con");
  auto final_state = std::make_shared<Matter>(true_pot, params);
  final_state->con2matter("product.con");
  auto init_path = helper_functions::neb_paths::linearPath(
      *initial, *final_state, params->nebImages);
  auto init_data = helper_functions::surrogate::getMidSlice(init_path);
  auto features = helper_functions::surrogate::get_features(init_data);
  auto targets = helper_functions::surrogate::get_targets(init_data, true_pot);
  // pybind11::scoped_interpreter guard{}; // Initialize the Python interpreter
  params->potential = PotType::GPR_Optim;
  auto surpot = helpers::create::makeSurrogatePotential(PotType::GPR_Optim,
                                                        params, initial);
  surpot->train_optimize(features, targets);
  auto matter = std::make_shared<Matter>(surpot, params);
  matter->con2matter("reactant.con");
  // auto [energy, forces, vari] = surpot->get_ef_var(
  //     init_path[0].getPositionsFree(), init_path[0].getAtomicNrsFree(),
  //     init_path[0].getCell());
  // auto execString = fmt::format("Got energy: {energy:}\n For
  // forces\n{forces:}\n With an energy variance of {variance:}\n",
  //                               fmt::arg("energy", energy),
  //                               fmt::arg("forces", fmt::streamed(forces)),
  //                               fmt::arg("variance", vari));
  // init_path[0].setPotential(surpot);
  auto execString = fmt::format(
      "Got energy: {energy:}\n For forces\n{forces:}\n With an energy variance "
      "of {variance:}\n",
      fmt::arg("energy", init_path[0].getPotentialEnergy()),
      fmt::arg("forces", fmt::streamed(init_path[0].getForcesFree())),
      fmt::arg("variance", init_path[0].getEnergyVariance()));
  std::cout << execString << "\n";
  return EXIT_SUCCESS;
}
