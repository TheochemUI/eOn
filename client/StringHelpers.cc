#include "StringHelpers.hpp"

#include <cassert>
#include <type_traits>

namespace helper_functions {
template <typename T>
std::vector<T> get_val_from_string(const std::string &line, std::optional<size_t> nelements) {
    assert(not line.empty());
    std::istringstream ss;
    std::vector<T> retval;
    auto elements{get_split_strings(line)};
    if (nelements.has_value()) {
        // Used to truncate if the number of elements is given
        assert(nelements > 0);
        elements.resize(nelements.value());
    }
    // If it is unsigned then use long double else T
    for (typename std::conditional<std::is_unsigned<T>::value, long double, T>::type tmp; auto elem : elements) {
        ss.str(elem); // elem is a string
        ss >> tmp;
        if (std::is_unsigned<T>::value) {
            assert(tmp > 0);
        }
        retval.push_back(tmp);
        ss.clear();
    }
    return retval;
}
// Instantiate explicitly
// This is useful since we don't want to support other types either
template std::vector<size_t> get_val_from_string(const std::string &, std::optional<size_t>);
template std::vector<double> get_val_from_string(const std::string &, std::optional<size_t>);

std::vector<std::string> get_split_strings(const std::string &line) {
    std::istringstream ss{line};
    std::vector<std::string> split_strings{std::istream_iterator<std::string>{ss},
                                           std::istream_iterator<std::string>()};
    return split_strings;
}
} // namespace helper_functions
