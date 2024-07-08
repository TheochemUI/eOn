/*
** This file is part of eON.
**
** SPDX-License-Identifier: BSD-3-Clause
**
** Copyright (c) 2010--present, eON Development Team
** All rights reserved.
**
** Repo:
** https://github.com/TheochemUI/eON
*/
#include "Matter.h"
#include "Parameters.h"
#include <cstdlib>

int main(void) {
  std::string confile("pos.con");
  auto parameters = std::make_shared<Parameters>();
  auto pot = helper_functions::makePotential(parameters);
  Matter *matter = new Matter(pot, parameters);
  matter->con2matter(confile);
  matter->writeTibble("rSysdat.txt"s);
  delete matter;
  return EXIT_SUCCESS;
}
