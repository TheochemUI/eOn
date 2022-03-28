/*
 * MorsePotTest.cpp
 *
 *  Created on: 27 Mar 2022
 *      Author: Rohit Goswami
 *     Company: University of Iceland
 */

#include <algorithm>

#include "../Potential.h"
#include "../potentials/Morse/Morse.h"

#include "../HelperFunctions.h"
#include "../GPRHelpers.h"
#include "../Matter.h"
#include "../Parameters.h"
#include "MorsePotTest.h"

namespace tests {

MorsePotTest::MorsePotTest() {
  // TODO Auto-generated constructor stub
}

MorsePotTest::~MorsePotTest() {
  // TODO Auto-generated destructor stub
}

TEST_F(MorsePotTest, TestMatter) {
  string confile("pos.con");
  Parameters *parameters = new Parameters;
  parameters->potential = "morse_pt";
  Matter *matter = new Matter(parameters);
  matter->con2matter(confile);
  Morse pot;
  matter->con2matter(confile);
  int nAtoms = matter->numberOfAtoms();
  auto posdata = matter->getPositions();
  auto celldat = matter->getCell();
  AtomMatrix forces = AtomMatrix::Constant(nAtoms, 3, 0);
  double *pos = posdata.data();
  double *frcs = forces.data();
  double *bx = celldat.data();
  double energy{0};
  pot.force(nAtoms, pos, nullptr, frcs, &energy, bx, 1);
  // TODO: Find a less hacky way
  AtomMatrix finForces{forces};
  for (int i = 0; i <nAtoms; i++){
    if(matter->getFixed(i)){
      finForces.row(i).setZero();
    }
  }
  EXPECT_EQ(finForces, matter->getForces())
      << "Forces do not match";
  EXPECT_EQ(energy, matter->getPotentialEnergy())
      << "Potential energy does not match";
  delete matter;
  delete parameters;
}

} /* namespace tests */
