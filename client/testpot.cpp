#include "Matter.h"
#include "Parameters.h"
#include "Potential.h"
#include <cstdlib>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <memory>
#include <pybind11/embed.h>

using namespace std::string_literals; // For ""s

int main(void) {
  Parameters *parameters = new Parameters;
  // parameters->potential = PotType::PYSURROGATE;
  // Potential *pot = helper_functions::makePotential(parameters);
  string confile("pos.con");
  auto params = std::make_shared<Parameters>();
  pybind11::scoped_interpreter guard{}; // Initialize the Python interpreter
  params->potential = PotType::PYSURROGATE;
  auto pot = helper_functions::makePotential(params);
  auto matter = std::make_unique<Matter>(pot, params);
  matter->con2matter(confile);
  auto [energy, forces, vari] = pot->get_ef(matter->getPositions(), matter->getAtomicNrs(), matter->getCell());
  auto execString =
    fmt::format("Got {energy:}\n{forces:}",
                  fmt::arg("energy", energy),
                  fmt::arg("forces", fmt::streamed(forces)));
  std::cout << execString << "\n";
  return EXIT_SUCCESS;
}
