#ifndef NEBJOBTEST_H_
#define NEBJOBGPRTEST_H_

#include <gtest/gtest.h>

#include <algorithm>

#include "../GPRHelpers.h"
#include "../HelperFunctions.h"
#include "../Log.h"
#include "../Matter.h"
#include "../NudgedElasticBandJob.h"
#include "../Parameters.h"
#include "../Potential.h"
#include "../potentials/Morse/Morse.h"

namespace tests {
    class NEBJobGPRTest : public ::testing::Test {
        public:
            NEBJobGPRTest();
            virtual ~NEBJobGPRTest();

    const long double threshold{1e-6};

    std::string reactantFilename;
    std::string productFilename;

    std::unique_ptr<Parameters> parameters;
    std::unique_ptr<Matter> initmatter;
    std::unique_ptr<Matter> finalmatter;
    std::unique_ptr<gpr::GaussianProcessRegression> gprfunc;
    std::vector<Matter> imgArray;
    gpr::Observation obspath;
    std::pair<gpr::AtomsConfiguration, gpr::Coord> config_data;
    AtomMatrix init_frcsref;
    double init_eref;
    std::function<bool(const gpr::EigenMatrix &lhs, const gpr::EigenMatrix &rhs)> comparer;
    gpr::AtomsConfiguration atoms_conf;
    };
} /* namespace tests */

#endif // NEBJOBGPRTEST_H_
