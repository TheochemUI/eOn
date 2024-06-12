#ifndef ENVHELPERS_H_
#define ENVHELPERS_H_

#include <string>

namespace helper_functions {
std::string get_value_from_env_or_param(const char *env_variable,
                                        const std::string &param_value,
                                        const std::string &default_value = "",
                                        const std::string &warning_message = "",
                                        const bool is_mandatory = false);
}

#endif // ENVHELPERS_H_
