#ifndef MATTERTEST_H
#define MATTERTEST_H

#include <gtest/gtest.h>
#include <string>

#include "../MatrixHelpers.hpp"
#include "../Matter.h"
#include "../Parameters.h"

using namespace std::string_literals; // For ""s

namespace tests {
class MatterTest : public ::testing::Test {
protected:
  double threshold;

public:
  MatterTest();
  virtual ~MatterTest();
};
} /* namespace tests */

#endif /* MATTERTEST_H */
