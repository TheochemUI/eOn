#include "BaseStructures.h"
#include "GPSurrogateJob.h"
#include "Matter.h"
#include "NudgedElasticBand.h"
#include "Parameters.h"
#include "Potential.h"
#include "SurrogatePotential.h"
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
  auto pot = helper_functions::makePotential(PotType::GPR_Optim, params);
  std::shared_ptr<SurrogatePotential> gp_pot =
      std::dynamic_pointer_cast<SurrogatePotential>(pot);

  if (!dynamic_cast<SurrogatePotential *>(gp_pot.get())) {
    throw "Ouch";
  }
  gp_pot->prepare(initial);
  gp_pot->train_optimize(features, targets);
  auto matter = std::make_shared<Matter>(gp_pot, params);
  matter->con2matter("reactant.con");
  auto [energy, forces, vari] = gp_pot->get_ef_var(
      init_path[2].getPositionsFree(), init_path[2].getAtomicNrsFree(),
      init_path[2].getCell());
  auto execString = fmt::format("Got {energy:}\n{forces:}\n{variance:}\n",
                                fmt::arg("energy", energy),
                                fmt::arg("forces", fmt::streamed(forces)),
                                fmt::arg("variance", fmt::streamed(vari)));
  std::cout << execString << "\n";
  return EXIT_SUCCESS;
}
