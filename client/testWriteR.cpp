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
#include <cstdlib>

using namespace std;

int main(void) {
  string confile("pos.con");
  Parameters parameters;
  auto pot = eonc::helpers::makePotential(parameters);
  Matter matter(pot, parameters);
  matter.con2matter(confile);
  matter.writeTibble("rSysdat.txt"s);
  return EXIT_SUCCESS;
}
