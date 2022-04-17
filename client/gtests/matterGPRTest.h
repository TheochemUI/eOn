#ifndef MATTERGPRTEST_H
#define MATTERGPRTEST_H

#include <gtest/gtest.h>

namespace tests {
    class matterGPRTest : public ::testing::Test {
        public:
            matterGPRTest();
            virtual ~matterGPRTest();
    string reactantFilename;
    string productFilename;
    const long double threshold{1e-3};

    Parameters params;
    };
} /* namespace tests */


#endif /* MATTERGPRTEST_H */
