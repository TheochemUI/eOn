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

bool VectorEquality(const VectorXd &lhs, const VectorXd &rhs) {
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

} /* namespace tests */
