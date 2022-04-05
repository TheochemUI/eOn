#ifndef NEBTEST_H_
#define NEBTEST_H_

#include <gtest/gtest.h>

namespace tests {
    class NEBTest : public ::testing::Test {
        public:
            NEBTest();
            virtual ~NEBTest();

    string reactantFilename;
    string productFilename;

    std::unique_ptr<Parameters> parameters;
    std::unique_ptr<Matter> initmatter;
    std::unique_ptr<Matter> finalmatter;
    };
} /* namespace tests */

#endif // NEBTEST_H_
