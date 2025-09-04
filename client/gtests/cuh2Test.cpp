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
/*
 * CuH2Test.cpp
 *
 *  Created on: 15 Nov 2022
 *      Author: Rohit Goswami
 *     Company: University of Iceland
 */

#include <algorithm>

#include "../Matter.h"
#include "../Parameters.h"
#include "CuH2Test.h"

namespace tests {

CuH2Test::CuH2Test() {
  // TODO Auto-generated constructor stub
}

CuH2Test::~CuH2Test() {
  // TODO Auto-generated destructor stub
}

TEST_F(CuH2Test, TestMatter) {
  string confile("pos.con");
  Parameters *parameters = new Parameters;
  parameters->potential = "cuh2_pot";
  Matter *matter = new Matter(parameters);
  matter->con2matter(confile);
  std::cout << matter->getPotentialEnergy();
  delete matter;
  delete parameters;
}

} /* namespace tests */
