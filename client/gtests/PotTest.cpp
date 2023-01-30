/*
 * Test.cpp
 *
 *  Created on: 30 Jan 2023
 *      Author: Rohit Goswami
 *     Company: University of Iceland
 */

#include "PotTest.h"
#include "../potentials/Morse/Morse.cpp"

#include <algorithm>

namespace tests {

PotTest::PotTest() {
    // TODO Auto-generated constructor stub
}

PotTest::~PotTest() {
    // TODO Auto-generated destructor stub
}

// Kanged from https://stackoverflow.com/a/39238772/1895378
bool MatrixEquality(const MatrixXd &lhs, const MatrixXd &rhs) {
  return lhs.isApprox(rhs, 1e-4);
}

bool VectorEquality(const VectorXd &lhs, const VectorXd &rhs) {
  return lhs.isApprox(rhs, 1e-4);
}

    TEST(PotTest, getName){
        string confile("pos.con");
        Parameters *parameters = new Parameters;
        parameters->potential = "lj";
        Potential* pot = Potential::getPotential(parameters);
        ASSERT_EQ(pot->getName(), "lj");
        parameters->potential = "morse_pt";
        Potential* pot2 = Potential::getPotential(parameters);
        ASSERT_EQ(pot2->getName(), "morse_pt");
    }

} /* namespace tests */
