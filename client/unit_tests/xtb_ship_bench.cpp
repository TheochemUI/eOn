// Time eOn linked XTBPot (ship) on water GFN2 — same geometry as rgpot bench.
#include "eon/Parameters.h"
#include "eon/Potential.h"
#include "eon/potentials/XTBPot/XTBPot.h"

#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

using SteadyClock = std::chrono::steady_clock;
using eonc::Parameters;
using eonc::PotType;

static const double water_pos[] = {0.00000000, 0.00000000,  0.11779000,
                                   0.00000000, 0.75545000,  -0.47116000,
                                   0.00000000, -0.75545000, -0.47116000};
static const int water_atmnrs[] = {8, 1, 1};
static const double water_box[9] = {100, 0, 0, 0, 100, 0, 0, 0, 100};

int main(int argc, char **argv) {
  std::string json_out;
  int warmup = 5;
  int iters = 50;
  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    if (a == "--json" && i + 1 < argc)
      json_out = argv[++i];
    else if (a == "--warmup" && i + 1 < argc)
      warmup = std::atoi(argv[++i]);
    else if (a == "--iters" && i + 1 < argc)
      iters = std::atoi(argv[++i]);
  }

  Parameters params;
  params.potential_options.potential = PotType::XTB;
  params.xtb_options.paramset = "GFN2xTB";
  params.xtb_options.acc = 1.0;
  params.xtb_options.elec_temperature = 300.0;
  params.xtb_options.maxiter = 250;
  params.xtb_options.charge = 0.0;
  params.xtb_options.uhf = 0;

  auto pot = std::make_shared<XTBPot>(params);
  double energy = 0;
  std::vector<double> forces(9, 0.0);
  auto once = [&] {
    pot->force(3, water_pos, water_atmnrs, forces.data(), &energy, nullptr,
               water_box);
  };

  for (int i = 0; i < warmup; ++i)
    once();
  double sum_ms = 0;
  for (int i = 0; i < iters; ++i) {
    auto t0 = SteadyClock::now();
    once();
    auto t1 = SteadyClock::now();
    sum_ms += std::chrono::duration<double, std::milli>(t1 - t0).count();
  }
  const double mean_ms = sum_ms / static_cast<double>(iters);
  std::cout << "eon_xtb_ship method=GFN2 water n_atoms=3 mean_ms=" << mean_ms
            << " energy_eV=" << energy << "\n";
  if (!json_out.empty()) {
    std::ofstream o(json_out);
    o << "{\n"
      << "  \"backend\": \"eon_linked_ship\",\n"
      << "  \"system\": \"water\",\n"
      << "  \"method\": \"GFN2xTB\",\n"
      << "  \"n_atoms\": 3,\n"
      << "  \"warmup\": " << warmup << ",\n"
      << "  \"iters\": " << iters << ",\n"
      << "  \"eon_linked_mean_ms\": " << mean_ms << ",\n"
      << "  \"energy_eV\": " << energy << "\n"
      << "}\n";
  }
  return 0;
}
