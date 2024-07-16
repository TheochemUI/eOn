#pragma once

#include "LJ/LJ.h"

namespace eonc::pot {
void from_toml(LJ::Params &p_a, const toml::node_view<const toml::node> &tbl_a);
}
