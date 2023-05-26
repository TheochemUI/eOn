#ifndef NEBTEST_H
#define NEBTEST_H

#include <gtest/gtest.h>
#include <string>

#include "../MatrixHelpers.hpp"
#include "../Matter.h"
#include "../Parameters.h"

using namespace std::string_literals; // For ""s

namespace tests {
class NEBTest : public ::testing::Test {
protected:
  Parameters *params;
  Matter *m1;
  Matter *m2;
  double threshold;
  void SetUp() override;
  void TearDown() override;

public:
  NEBTest();
  virtual ~NEBTest();
};
} /* namespace tests */

#endif /* NEBTEST_H */
