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
  Parameters *params;
  Matter *m1;
  double threshold;
  void SetUp() override;
  void TearDown() override;

public:
  MatterTest();
  virtual ~MatterTest();
};
} /* namespace tests */

#endif /* MATTERTEST_H */
