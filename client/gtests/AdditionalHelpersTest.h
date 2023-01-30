#ifndef ADDITIONALHELPERSTEST_H
#define ADDITIONALHELPERSTEST_H

#include <gtest/gtest.h>
#include <string>

namespace tests {
class AdditionalHelpersTest : public ::testing::Test {
public:
    AdditionalHelpersTest();
    virtual ~AdditionalHelpersTest();
    std::string number_string{"1 2.6 4 5"s};
};
} /* namespace tests */

#endif /* ADDITIONALHELPERSTEST_H */
