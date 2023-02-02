/*
 * Test.cpp
 *
 *  Created on: 23 Feb 2022
 *      Author: Rohit Goswami
 *     Company: University of Iceland
 */

#include "NEBTest.h"
#include "../potentials/Morse/Morse.cpp"
#include "NudgedElasticBand.h"
#include "Log.h"

#include <algorithm>

using namespace std::placeholders;

namespace tests {

NEBTest::NEBTest()
    : params{new Parameters}, m1{new Matter(params)}, m2{new Matter(params)},
      threshold{1e-6} {}

NEBTest::~NEBTest() {
  delete params;
  delete m1;
}

void NEBTest::SetUp() {
  std::string confile("pos.con");
  m1->con2matter(confile);
  m2->con2matter(confile);
  m2->setPositions(m1->getPositions() * 3);
  log_init(params, "blah.log");
}

void NEBTest::TearDown() {}

TEST_F(NEBTest, TestCreation) {
  NudgedElasticBand* neb = new NudgedElasticBand(m1, m2, params);
  int status = neb->compute();
  ASSERT_EQ(status, 3);
}

// TEST_F(NEBTest, SetGetAtomicNrs) {
//   VectorXi _atmnrs{{8, 8, 6, 6, 6, 6, 1, 1, 1, 1, 1, 1, 16}};
//   VectorXi _atmnrs2{{16, 16, 12, 12, 12, 12, 2, 2, 2, 2, 2, 2, 32}};
//   auto vecEq =
//       std::bind(helper_functions::eigenEquality<VectorXi>, _1, _2, threshold);
//   ASSERT_PRED2(vecEq, _atmnrs, m1->getAtomicNrs());
//   for (auto &atmnr : _atmnrs) {
//     atmnr *= 2;
//   }
//   m1->setAtomicNrs(_atmnrs);
//   ASSERT_PRED2(vecEq, _atmnrs2, m1->getAtomicNrs());
// }

// TEST_F(NEBTest, SetPotential) {
//   params->potential = "lj";
//   double m1_ipot = m1->getPotentialEnergy();
//   params->potential = "morse_pt";
//   Potential *pot{Potential::getPotential(params)};
//   ASSERT_NE(m1->getPotential(), pot);
//   m1->setPotential(pot);
//   ASSERT_EQ(pot->getName(), "morse_pt");
//   ASSERT_EQ(m1->getPotential(), pot);
//   double m1_fpot = m1->getPotentialEnergy();
//   ASSERT_NE(m1_ipot, m1_fpot);
// }

} /* namespace tests */
