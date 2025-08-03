#pragma once

#ifdef EON_CHECKS
#include "fmt/core.h"
#endif

#include "magic_enum/magic_enum.hpp"
#include "thirdparty/toml.hpp"
#include <string>

namespace eonc {
// TODO(rg): Maybe add a layer of indirection and allow other kinds of inputs...
toml::table loadTOML(const std::string &fname);
void config_section(const toml::table &conf, const std::string_view key);

template <typename T>
std::optional<T> get_enum_toml(const toml::node_view<const toml::node> &conf) {
#ifdef EON_CHECKS
  if (!conf.is_string()) {
    throw std::logic_error("Node not found in TOML");
  }
#endif

  auto cval = conf.value<std::string>();
  auto enum_value =
      magic_enum::enum_cast<T>(cval.value(), magic_enum::case_insensitive);
#ifdef EON_CHECKS
  if (!enum_value.has_value()) {
    throw std::runtime_error(
        fmt::format("{} does not map to a valid option", cval.value()));
  }
#endif

  return enum_value;
}
} // namespace eonc
