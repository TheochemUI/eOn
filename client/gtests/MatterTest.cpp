/*
 * Test.cpp
 *
 *  Created on: 23 Feb 2022
 *      Author: Rohit Goswami
 *     Company: University of Iceland
 */

#include "MatterTest.h"

#include <algorithm>

namespace tests {

MatterTest::MatterTest() {
    // TODO Auto-generated constructor stub
}

MatterTest::~MatterTest() {
    // TODO Auto-generated destructor stub
}

TEST_F(MatterTest, TestCell) {
    Matter *m1 = new Matter(&p);
    m1->con2matter(fname);
    long _nAtoms{2+4+6+1};
    AtomMatrix _positions;
    AtomMatrix _velocities;
    AtomMatrix _forces;
    AtomMatrix _biasForces;
    // VectorXd _masses{{15.99, 12.011, 1.008, 32.065}}; 3.4.0
    VectorXd _masses(4);
    VectorXi _atomicNrs;
    VectorXi _isFixed;
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
    _masses << 15.99, 12.011, 1.008, 32.065;
    // clang-format on
    ASSERT_EQ(m1->numberOfAtoms(), _nAtoms);
    // ASSERT_EQ(m1->getPositions(), _positions);
    // ASSERT_EQ(m1->getPositionsV(), VectorXd::Map(_positions.data(), 3 * _nAtoms));
    // ASSERT_EQ(m1->getCell(), _cell);
    // ASSERT_EQ(m1->getMasses(), _masses);
    delete m1;
}

} /* namespace tests */
