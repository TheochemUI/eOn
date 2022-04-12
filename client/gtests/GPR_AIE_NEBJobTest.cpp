/*
 * GPR_AIE_NEBJobTest.cpp
 *
 *  Created on: 6 April 2022
 *      Author: Rohit Goswami
 *     Company: University of Iceland
 */

#include <algorithm>

#include "../Potential.h"
#include "../potentials/Morse/Morse.h"
#include "../GPR_AIE_NEBJob.h"
#include "../Job.h"

#include "../HelperFunctions.h"
#include "../GPRHelpers.h"
#include "../Log.h"
#include "../Matter.h"
#include "../Parameters.h"
#include "GPR_AIE_NEBJobTest.h"

namespace tests {

GPR_AIE_NEBJobTest::GPR_AIE_NEBJobTest() {
    parameters = std::make_unique<Parameters>();
    parameters->potential = "morse_pt";
    parameters->nebImages = 7;
    parameters->LogPotential = false;
    log_init(parameters.get(), (char *)"test.log");
}

GPR_AIE_NEBJobTest::~GPR_AIE_NEBJobTest() {
}

TEST_F(GPR_AIE_NEBJobTest, TestRun) {
        auto job = new GPR_AIE_NEBJob(parameters.get());
        job->run();
        delete job;
}

} /* namespace tests */
