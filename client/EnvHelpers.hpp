#ifndef ENVHELPERS_H_
#define ENVHELPERS_H_

#include <string>

namespace helper_functions {
std::string get_value_from_param_or_env(const std::string &param_value,
                                        const char *env_variable,
                                        const std::string &default_value = "",
                                        const std::string &warning_message = "",
                                        const bool is_mandatory = false);
}

#endif // ENVHELPERS_H_
