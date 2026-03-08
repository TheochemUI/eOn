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
#include "PyGuard.h"
#include <cstdlib>
#include <iostream>
#include <memory>
#include <sstream>

using namespace std;
using namespace std::string_literals; // For ""s

int main(void) {
  string confile("pos.con");
  auto params = std::make_shared<Parameters>();
  eonc::ensure_interpreter();
  params->potential_options.potential = PotType::CatLearn;
  auto pot = eonc::helpers::makePotential(params);
  auto matter = std::make_unique<Matter>(pot, params);
  matter->con2matter(confile);
  auto [energy, forces] = pot->get_ef(
      matter->getPositions(), matter->getAtomicNrs(), matter->getCell());
  std::ostringstream oss;
  oss << "Got " << energy << "\n" << forces << "\n";
  std::cout << oss.str();
  return EXIT_SUCCESS;
}
