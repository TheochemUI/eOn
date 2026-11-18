/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
*/
#include "eon/ConFileIO.h"
#include "eon/Matter.h"
#include "eon/Parameters.h"
#include "eon/Potential.h"
#include "eon/potentials/LAMMPS/LAMMPSPot.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

int main(int argc, char **argv) {
  if (argc != 4) {
    std::cerr << "usage: cu_force_benchmark MODE POS_CON REPEATS\n";
    return 2;
  }
  const std::string mode = argv[1];
  const std::string con = argv[2];
  const int repeats = std::stoi(argv[3]);
  if (repeats < 6) {
    throw std::invalid_argument("REPEATS must be at least 6");
  }

  eonc::Parameters params;
  eonc::ParametersLoadAccess::potential_options(params).potential =
      eonc::PotType::LAMMPS;
  std::shared_ptr<eonc::Potential> potential;
  if (mode == "worker" || mode == "compare") {
    potential = std::make_shared<LAMMPSPot>(
        params, LAMMPSPot::Execution::Worker);
  } else if (mode == "library") {
    potential = std::make_shared<LAMMPSPot>(params);
  } else if (mode == "rgpot") {
    potential =
        eonc::helpers::sharePotential(eonc::helpers::makePotential(params));
  } else {
    throw std::invalid_argument(
        "MODE must be worker, library, rgpot, or compare");
  }

  eonc::Matter matter(potential, params);
  if (!eonc::io::io_ok(matter.con2matter(con))) {
    throw std::runtime_error("could not read " + con);
  }
  const long n = matter.numberOfAtoms();
  if (n != 4000) {
    throw std::runtime_error("benchmark input must contain exactly 4000 atoms");
  }
  const AtomMatrix positions = matter.getPositionsCopy();
  const VectorXi atomicNumbers = matter.getAtomicNrs();
  const Matrix3d box = matter.getCell();
  AtomMatrix forces = AtomMatrix::Zero(n, 3);
  double energy = 0.0;
  double variance = 0.0;

  if (mode == "compare") {
    auto native =
        eonc::helpers::sharePotential(eonc::helpers::makePotential(params));
    AtomMatrix nativeForces = AtomMatrix::Zero(n, 3);
    double nativeEnergy = 0.0;
    double nativeVariance = 0.0;
    potential->force(n, positions.data(), atomicNumbers.data(), forces.data(),
                     &energy, &variance, box.data());
    native->force(n, positions.data(), atomicNumbers.data(),
                  nativeForces.data(), &nativeEnergy, &nativeVariance,
                  box.data());
    const double energyDifference = std::abs(nativeEnergy - energy);
    const double forceDifference =
        (nativeForces - forces).cwiseAbs().maxCoeff();
    std::cout << std::setprecision(12) << "energy_difference_eV,"
              << energyDifference << "\nforce_linf_difference_eV_per_A,"
              << forceDifference << '\n';
    if (energyDifference > 1e-8 || forceDifference > 1e-8) {
      return 1;
    }
    return 0;
  }

  for (int warmup = 0; warmup < 2; ++warmup) {
    potential->force(n, positions.data(), atomicNumbers.data(), forces.data(),
                     &energy, &variance, box.data());
  }

  std::cout << "mode,repeat,atoms,latency_ms,energy_eV,force_linf_eV_per_A\n";
  for (int repeat = 0; repeat < repeats; ++repeat) {
    const auto start = std::chrono::steady_clock::now();
    potential->force(n, positions.data(), atomicNumbers.data(), forces.data(),
                     &energy, &variance, box.data());
    const double milliseconds =
        std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - start)
            .count();
    const double forceMax = forces.cwiseAbs().maxCoeff();
    if (!std::isfinite(energy) || energy == 0.0 ||
        !std::isfinite(forceMax)) {
      throw std::runtime_error("non-finite or zero benchmark result");
    }
    std::cout << mode << ',' << repeat << ',' << n << ',' << std::setprecision(12)
              << milliseconds << ',' << energy << ',' << forceMax << '\n';
  }
}
