#ifndef GPR_AIE_NEBJOBTEST_H_
#define GPR_AIE_NEBJOBTEST_H_

#include <gtest/gtest.h>

namespace tests {
    class GPR_AIE_NEBJobTest : public ::testing::Test {
        public:
            GPR_AIE_NEBJobTest();
            virtual ~GPR_AIE_NEBJobTest();
    };
    std::unique_ptr<Parameters> parameters;
} /* namespace tests */

#endif // GPR_AIE_NEBJOBTEST_H_
