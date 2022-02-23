/*
 * StringHelpersTest.cpp
 *
 *  Created on: 23 Feb 2022
 *      Author: Rohit Goswami
 *     Company: University of Iceland
 */

#include "../StringHelpers.hpp"

#include "StringHelpersTest.h"
#include <vector>

using namespace std::string_literals; // For ""s

namespace tests {

StringHelpersTest::StringHelpersTest() {
    // TODO Auto-generated constructor stub
}

StringHelpersTest::~StringHelpersTest() {
    // TODO Auto-generated destructor stub
}

TEST_F(StringHelpersTest, TestSplitStrings) {
    auto split_strings = helper_functions::get_split_strings(number_string);
    ASSERT_EQ(split_strings.size(), 4) << "Size mismatch after split"s;
    EXPECT_EQ(split_strings[0], "1"s);
    EXPECT_EQ(split_strings[1], "2.6"s);
    EXPECT_EQ(split_strings[2], "4"s);
    EXPECT_EQ(split_strings[3], "5"s);
}

TEST_F(StringHelpersTest, TestValsBasic) {
    // Convert to doubles
    auto split_strings = helper_functions::get_val_from_string<double>(number_string);
    ASSERT_EQ(split_strings.size(), 4) << "Size mismatch after split"s;
    EXPECT_EQ(split_strings[0], 1.0);
    EXPECT_EQ(split_strings[1], 2.60);
    EXPECT_EQ(split_strings[2], 4.0);
    EXPECT_EQ(split_strings[3], 5.0);
    // Test with size_t
    auto split_strings_size_t = helper_functions::get_val_from_string<size_t>(number_string);
    ASSERT_EQ(split_strings_size_t.size(), 4) << "Size mismatch after split"s;
    EXPECT_EQ(split_strings_size_t[0], 1);
    EXPECT_EQ(split_strings_size_t[1], 2);
    EXPECT_EQ(split_strings_size_t[2], 4);
    EXPECT_EQ(split_strings_size_t[3], 5);
}

TEST_F(StringHelpersTest, TestValsTrunc) {
    // Convert to doubles
    auto split_strings = helper_functions::get_val_from_string<double>(number_string, 2);
    ASSERT_EQ(split_strings.size(), 2);
    EXPECT_EQ(split_strings[0], 1.0);
    EXPECT_EQ(split_strings[1], 2.60);
    // Test with size_t
    auto split_strings_size_t = helper_functions::get_val_from_string<size_t>(number_string, 3);
    ASSERT_EQ(split_strings_size_t.size(), 3);
    EXPECT_EQ(split_strings_size_t[0], 1);
    EXPECT_EQ(split_strings_size_t[1], 2);
    EXPECT_EQ(split_strings_size_t[2], 4);
}

TEST_F(StringHelpersTest, TestEmpty) {
    std::string number_string{""s};
    auto split_strings = helper_functions::get_split_strings(number_string);
    ASSERT_EQ(split_strings.size(), 0);
    EXPECT_DEATH(helper_functions::get_val_from_string<double>(number_string, 1),
                 "(not line.empty())"s);
    EXPECT_DEATH(helper_functions::get_val_from_string<size_t>(number_string),
                 "(not line.empty())"s);
}

TEST_F(StringHelpersTest, TestNeg) {
    std::string tester{"1 -1"s};
    auto split_strings = helper_functions::get_split_strings(tester);
    ASSERT_EQ(split_strings.size(), 2);
    // This needs to fail, since size_t is unsigned
    EXPECT_DEATH(helper_functions::get_val_from_string<size_t>(tester),
                 "Can't represent negative numbers with an unsigned type"s);
    // Should be OK, due to truncation
    EXPECT_EQ(helper_functions::get_val_from_string<size_t>(tester, 1)[0], 1);
    // Fine otherwise
    std::vector<double> exp_vals = {1.0, -1.0};
    EXPECT_EQ(helper_functions::get_val_from_string<double>(tester), exp_vals);
}

} /* namespace tests */
