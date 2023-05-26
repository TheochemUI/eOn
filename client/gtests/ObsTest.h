#ifndef OBSTEST_H
#define OBSTEST_H

#include <gtest/gtest.h>

namespace tests {
class ObsTest : public ::testing::Test {
public:
  ObsTest();
  virtual ~ObsTest();
};
} /* namespace tests */

#endif /* OBSTEST_H */
