/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
**
** Copyright (c) 2010--present, eOn Development Team
** All rights reserved.
**
** Repo:
** https://github.com/TheochemUI/eOn
*/
#include "Matter.h"
#include "Parameters.h"
#include "Potential.h"
#include <cstdlib>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <iostream>
#include <memory>
#include <pybind11/embed.h>

using namespace std::string_literals; // For ""s

int main(void) {
  std::string confile("pos.con");
  auto params = std::make_shared<Parameters>();
  pybind11::scoped_interpreter guard{}; // Initialize the Python interpreter
  params->pot.potential = PotType::CatLearn;
  auto pot = helper_functions::makePotential(params);
  auto matter = std::make_unique<Matter>(pot, params);
  matter->con2matter(confile);
  auto [energy, forces] = pot->get_ef(
      matter->getPositions(), matter->getAtomicNrs(), matter->getCell());
  auto execString =
      fmt::format("Got {energy:}\n{forces:}", fmt::arg("energy", energy),
                  fmt::arg("forces", fmt::streamed(forces)));
  std::cout << execString << "\n";
  return EXIT_SUCCESS;
}
