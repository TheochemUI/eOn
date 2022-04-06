#ifndef GPRPOTTEST_H
#define GPRPOTTEST_H

#include <gtest/gtest.h>

namespace tests {
    class GPRPotTest : public ::testing::Test {
        public:
            GPRPotTest();
            virtual ~GPRPotTest();
    const long double threshold{1e-6};
    string reactantFilename;
    string productFilename;

    std::unique_ptr<Parameters> parameters;
    std::unique_ptr<Parameters> gprparameon;
    std::unique_ptr<Matter> initmatter;
    std::unique_ptr<Matter> finalmatter;
    std::unique_ptr<gpr::GaussianProcessRegression> gprfunc;
    };
} /* namespace tests */


#endif /* GPRPOTTEST_H */
