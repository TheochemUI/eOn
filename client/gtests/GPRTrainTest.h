#ifndef GPRTRAINTEST_H
#define GPRTRAINTEST_H

#include <gtest/gtest.h>

namespace tests {
    class GPRTrainTest : public ::testing::Test {
        public:
            GPRTrainTest();
            virtual ~GPRTrainTest();
    long double threshold;
    string reactantFilename;
    string productFilename;

    std::unique_ptr<Parameters> parameters;
    std::unique_ptr<Matter> initmatter;
    std::unique_ptr<Matter> finalmatter;
    };
} /* namespace tests */


#endif /* GPRTRAINTEST_H */
