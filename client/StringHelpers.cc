#include "StringHelpers.hpp"

#include <cassert>
#include <type_traits>

namespace helper_functions {
template <typename T>
std::vector<T> get_val_from_string(const std::string &line,
                                   std::optional<size_t> nelements) {
  assert(not line.empty());
  std::vector<T> retval;
  const bool b_isunsigned{std::is_unsigned<T>::value};
  auto elements{get_split_strings(line)};
  if (nelements.has_value()) {
    // Used to truncate if the number of elements is given
    assert(nelements > 0);
    elements.resize(nelements.value());
  }
  // If it is unsigned then use long double else T
  for (typename std::conditional<b_isunsigned, long double, T>::type tmp;
       auto elem : elements) {
    if (not isNumber(elem)) {
      continue;
    }
    std::istringstream ss{
        elem}; // instead of {ss.str(elem); ss >> tmp; ss.clear();}
    ss >> tmp;
    if (b_isunsigned and tmp < 0) {
      std::cerr
          << "Can't represent negative numbers with an unsigned type, bailing on "s
          << tmp << "\n";
      assert(tmp > 0);
    }
    retval.push_back(tmp);
  }
  return retval;
}
// Instantiate explicitly
// This is useful since we don't want to support other types either
template std::vector<size_t> get_val_from_string(const std::string &,
                                                 std::optional<size_t>);
template std::vector<double> get_val_from_string(const std::string &,
                                                 std::optional<size_t>);

std::vector<std::string> get_split_strings(const std::string &line) {
  std::istringstream ss{line};
  std::vector<std::string> split_strings{std::istream_iterator<std::string>{ss},
                                         std::istream_iterator<std::string>()};
  return split_strings;
}

bool isNumber(const std::string &token) {
  return std::regex_match(
      token, std::regex(("((\\+|-)?[[:digit:]]+)(\\.(([[:digit:]]+)?))?")));
}
} // namespace helper_functions
