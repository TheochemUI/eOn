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

    std::unique_ptr<GPRMatter> ptrGPRM;
    };
} /* namespace tests */


#endif /* MATTERGPRTEST_H */
