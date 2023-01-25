#ifndef MATTERTEST_H
#define MATTERTEST_H

#include <gtest/gtest.h>
#include <string>

#include "../Matter.h"
#include "../Parameters.h"

using namespace std::string_literals; // For ""s

namespace tests {
class MatterTest : public ::testing::Test {
public:
    MatterTest();
    virtual ~MatterTest();
    std::string fname {"pos.con"s};
    Parameters p;
};
} /* namespace tests */

#endif /* MATTERTEST_H */
