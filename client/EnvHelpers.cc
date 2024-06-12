#include "EnvHelpers.hpp"
#include <spdlog/spdlog.h>

namespace helper_functions {

std::string get_value_from_param_or_env(const std::string &param_value,
                                        const char *env_variable,
                                        const std::string &default_value,
                                        const std::string &warning_message,
                                        const bool is_mandatory) {
  if (!param_value.empty()) {
    return param_value;
  }

  const char *env_value = std::getenv(env_variable);
  if (env_value != nullptr) {
    return std::string(env_value);
  }

  if (is_mandatory) {
    throw std::runtime_error(
        "Environment variable " + std::string(env_variable) +
        " is not set and no parameter value provided. Please set it in the "
        "configuration or as an environment variable.\n");
  }

  if (!default_value.empty() && !warning_message.empty()) {
    SPDLOG_WARN(warning_message);
  }

  return default_value;
}

} // namespace helper_functions
