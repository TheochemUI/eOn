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

namespace eonc {



namespace helpers {
/**
 * \brief Parse a string into values
 *
 * @param line A thing to be parsed
 */
template <typename T>
std::vector<T>
get_val_from_string(std::string_view line,
                    std::optional<size_t> nelements = std::nullopt);
/**
 * \brief Split a string into constituent strings
 *
 * Based on https://www.fluentcpp.com/2017/04/21/how-to-split-a-string-in-c/
 *
 * @param line A thing to be parsed
 */
std::vector<std::string> get_split_strings(std::string_view line);
/**
 * \brief Figure out if a string has a number in it
 *
 * Based on
 *
 */
bool isNumber(std::string_view token);
} // namespace helpers

} // namespace eonc
