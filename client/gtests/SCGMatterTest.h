#ifndef SCGMATTERTEST_H
#define SCGMATTERTEST_H

#include <gtest/gtest.h>

namespace tests {
    class SCGMatterTest : public ::testing::Test {
        public:
            SCGMatterTest();
            virtual ~SCGMatterTest();

            funcmin::SCG scg;
            double threshold;  // Epsilon for comparison operators.
    };
} /* namespace tests */


#endif /* SCGMATTERTEST_H */
