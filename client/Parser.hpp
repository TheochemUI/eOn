#pragma once

#include <string>
#include "thirdparty/toml.hpp"

namespace eonc {
// TODO(rg): Maybe add a layer of indirection and allow other kinds of inputs...
toml::table loadTOML(const std::string fname);
} // namespace eonc
