#ifndef GPRPOTTEST_H
#define GPRPOTTEST_H

#include <algorithm>
#include <gtest/gtest.h>

#include "../Potential.h"
#include "../potentials/Morse/Morse.h"
#include "../potentials/GPRPotential/GPRPotential.h"

#include "../HelperFunctions.h"
#include "../GPRHelpers.h"
#include "../Matter.h"
#include "../Parameters.h"
#include "../subprojects/gprdimer/gpr/ml/GaussianProcessRegression.h"


namespace tests {
    class GPRPotTest : public ::testing::Test {
        public:
            GPRPotTest();
            virtual ~GPRPotTest();
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
    };
} /* namespace tests */


#endif /* GPRPOTTEST_H */
