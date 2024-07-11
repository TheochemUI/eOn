#pragma once

#include "thirdparty/toml.hpp"
#include <string>

namespace eonc {
// TODO(rg): Maybe add a layer of indirection and allow other kinds of inputs...
toml::table loadTOML(const std::string& fname);
void config_section(const toml::table &conf, const std::string_view key);
} // namespace eonc
