/*
 * Test.cpp
 *
 *  Created on: 23 Feb 2022
 *      Author: Rohit Goswami
 *     Company: University of Iceland
 */

#include "MatterTest.h"
#include "../potentials/Morse/Morse.cpp"

#include <algorithm>

using namespace std::placeholders;

namespace tests {

  MatterTest::MatterTest() : params{new Parameters}, m1{new Matter(params)}, threshold{1e-6} {}

MatterTest::~MatterTest() {
  delete params;
  delete m1;
}

void MatterTest::SetUp() {
  std::string confile("pos.con");
  m1->con2matter(confile);
}

void MatterTest::TearDown() {}

TEST_F(MatterTest, TestCell) {
  Matrix3d _cell;
  Matrix3d _cellInverse;
  auto matEq = std::bind(helper_functions::eigenEquality<Matrix3d>, _1, _2, threshold);
  // clang-format off
    _cell << // Comma initialized
        25.0, 0.0, 0.0,
        0.0, 25.0, 0.0,
        0.0, 0.0, 25.0;
  // clang-format on
  ASSERT_PRED2(matEq, _cell, m1->getCell());
}

TEST_F(MatterTest, SetGetAtomicNrs) {
  VectorXi _atmnrs {{ 8, 8, 6, 6, 6, 6, 1, 1, 1, 1, 1, 1, 16}};
  VectorXi _atmnrs2 {{ 16, 16, 12, 12, 12, 12, 2, 2, 2, 2, 2, 2, 32 }};
  auto vecEq = std::bind(helper_functions::eigenEquality<VectorXi>, _1, _2, threshold);
  ASSERT_PRED2(vecEq, _atmnrs, m1->getAtomicNrs());
  for (auto& atmnr : _atmnrs){
    atmnr *= 2;
  }
  m1->setAtomicNrs(_atmnrs);
  ASSERT_PRED2(vecEq, _atmnrs2, m1->getAtomicNrs());
}

} /* namespace tests */
