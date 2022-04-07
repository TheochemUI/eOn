#ifndef NEBJOBTEST_H_
#define NEBJOBGPRTEST_H_

#include <gtest/gtest.h>

namespace tests {
    class NEBJobGPRTest : public ::testing::Test {
        public:
            NEBJobGPRTest();
            virtual ~NEBJobGPRTest();

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

#endif // NEBJOBGPRTEST_H_
