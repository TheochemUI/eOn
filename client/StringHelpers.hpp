/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
**
** Copyright (c) 2010--present, eOn Development Team
** All rights reserved.
**
** Repo:
** https://github.com/TheochemUI/eOn
*/
#pragma once

#include <optional>
#include <regex>

#include "Matter.h"
#include "Parameters.h"

namespace helper_functions {
/**
 * \brief Parse a string into values
 *
 * @param std::string A thing to be parsed
 */
template <typename T>
std::vector<T>
get_val_from_string(const std::string &line,
                    std::optional<size_t> nelements = std::nullopt);
/**
 * \brief Split a string into constituent strings
 *
 * Based on https://www.fluentcpp.com/2017/04/21/how-to-split-a-string-in-c/
 *
 * @param std::string A thing to be parsed
 */
std::vector<std::string> get_split_strings(const std::string &line);
/**
 * \brief Figure out if a string has a number in it
 *
 * Based on
 *
 */
bool isNumber(const std::string &token);
} // namespace helper_functions
