#ifndef GPRTRAINTEST_H
#define GPRTRAINTEST_H

#include <gtest/gtest.h>
#include "../GPRHelpers.h"
#include "../HelperFunctions.h"
#include "../Matter.h"
#include "../Parameters.h"
#include "../Potential.h"
#include "../potentials/GPRPotential/GPRPotential.h"
#include "../potentials/Morse/Morse.h"
#include "../subprojects/gprdimer/gpr/ml/GaussianProcessRegression.h"

namespace tests {
    class GPRTrainTest : public ::testing::Test {
        public:
            GPRTrainTest();
            virtual ~GPRTrainTest();
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


#endif /* GPRTRAINTEST_H */
