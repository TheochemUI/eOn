#ifndef GPRTRAINTEST_H
#define GPRTRAINTEST_H

#include <gtest/gtest.h>

namespace tests {
    class GPRTrainTest : public ::testing::Test {
        public:
            GPRTrainTest();
            virtual ~GPRTrainTest();
    const long double threshold{1e-6};
    string reactantFilename;
    string productFilename;

    std::unique_ptr<Parameters> parameters;
    std::unique_ptr<Matter> initmatter;
    std::unique_ptr<Matter> finalmatter;
    std::unique_ptr<gpr::GaussianProcessRegression> gprfunc;
    };
} /* namespace tests */


#endif /* GPRTRAINTEST_H */
