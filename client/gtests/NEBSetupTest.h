#ifndef NEBSETUPTEST_H_
#define NEBSETUPTEST_H_

#include <gtest/gtest.h>

namespace tests {
    class NEBSetupTest : public ::testing::Test {
        public:
            NEBSetupTest();
            virtual ~NEBSetupTest();

    string reactantFilename;
    string productFilename;

    std::unique_ptr<Parameters> parameters;
    std::unique_ptr<Matter> initmatter;
    std::unique_ptr<Matter> finalmatter;
    };
} /* namespace tests */

#endif // NEBSETUPTEST_H_
