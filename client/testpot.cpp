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
  std::unique_ptr<Parameters> params = std::make_unique<Parameters>();
  pybind11::scoped_interpreter guard{}; // Initialize the Python interpreter
  params->potential = PotType::PYSURROGATE;
  Potential *pot = helper_functions::makePotential(params.get());
  std::unique_ptr<Matter> matter = std::make_unique<Matter>(params.get());
  matter->con2matter(confile);
  auto [energy, forces] = pot->get_ef(matter->getPositions(), matter->getAtomicNrs(), matter->getCell());
  auto execString =
    fmt::format("Got {energy:}\n{forces:}",
                  fmt::arg("energy", energy),
                  fmt::arg("forces", fmt::streamed(forces)));
  std::cout << execString << "\n";
  return EXIT_SUCCESS;
}
