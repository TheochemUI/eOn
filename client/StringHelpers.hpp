#ifndef STRINGHELPERS_H
#define STRINGHELPERS_H

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
#endif /* STRINGHELPERS_H */
