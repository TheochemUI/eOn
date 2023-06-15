#ifndef POTTEST_H
#define POTTEST_H

#include <gtest/gtest.h>
#include <string>

#include "../MatrixHelpers.hpp"
#include "../Matter.h"
#include "../Parameters.h"

using namespace std::string_literals; // For ""s

namespace tests {
class PotTest : public ::testing::Test {
protected:
  std::shared_ptr<Parameters> params;
  std::shared_ptr<Matter> m1;
  std::shared_ptr<Potential> pot_default;
  double threshold;
  void SetUp() override;
  void TearDown() override;

public:
  PotTest();
  virtual ~PotTest();
};
} /* namespace tests */

#endif /* POTTEST_H */
