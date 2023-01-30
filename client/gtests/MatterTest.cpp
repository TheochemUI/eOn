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

namespace tests {

MatterTest::MatterTest() : params{new Parameters}, m1{new Matter(params)} {}

MatterTest::~MatterTest() {
  delete params;
  delete m1;
}

void MatterTest::SetUp() {
  std::string confile("pos.con");
  m1->con2matter(confile);
}

void MatterTest::TearDown() {}

// Kanged from https://stackoverflow.com/a/39238772/1895378
bool MatrixEquality(const MatrixXd &lhs, const MatrixXd &rhs) {
  return lhs.isApprox(rhs, 1e-4);
}

bool VectorEquality(const VectorXi &lhs, const VectorXi &rhs) {
  return lhs.isApprox(rhs, 1e-4);
}

TEST_F(MatterTest, TestCell) {
  Matrix3d _cell;
  Matrix3d _cellInverse;
  // clang-format off
    _cell << // Comma initialized
        25.0, 0.0, 0.0,
        0.0, 25.0, 0.0,
        0.0, 0.0, 25.0;
    _cellInverse << // Comma initialized
        0.04, 0.0, 0.0,
        0.0, 0.04, 0.0,
        0.0, 0.0, 0.04;
  // clang-format on
  ASSERT_PRED2(MatrixEquality, _cell, m1->getCell());
  ASSERT_PRED2(MatrixEquality, _cellInverse, m1->getCell().inverse());
}

TEST_F(MatterTest, SetGetAtomicNrs) {
  VectorXi _atmnrs {{ 8, 8, 6, 6, 6, 6, 1, 1, 1, 1, 1, 1, 16}};
  VectorXi _atmnrs2 {{ 16, 16, 12, 12, 12, 12, 2, 2, 2, 2, 2, 2, 32 }};
  ASSERT_PRED2(VectorEquality, _atmnrs, m1->getAtomicNrs());
  for (auto& atmnr : _atmnrs){
    atmnr *= 2;
  }
  m1->setAtomicNrs(_atmnrs);
  ASSERT_PRED2(VectorEquality, _atmnrs2, m1->getAtomicNrs());
}

} /* namespace tests */
