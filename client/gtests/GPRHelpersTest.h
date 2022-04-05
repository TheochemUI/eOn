#ifndef GPRHELPERSTEST_H
#define GPRHELPERSTEST_H

#include <gtest/gtest.h>

namespace tests {
    class GPRHelpersTest : public ::testing::Test {
        public:
            GPRHelpersTest();
            virtual ~GPRHelpersTest();
    string reactantFilename;
    string productFilename;

    std::unique_ptr<Parameters> parameters;
    std::unique_ptr<Matter> initmatter;
    std::unique_ptr<Matter> finalmatter;
    };
} /* namespace tests */


#endif /* GPRHELPERSTEST_H */
