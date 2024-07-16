#include "ParseTOML.hpp"

namespace eonc::pot {

void from_toml(LJ::Params &p_a, const toml::node_view<const toml::node> &tbl) {
  p_a.u0 = tbl["u0"].value_or(p_a.u0);
  p_a.cutoff_R = tbl["cutoff"].value_or(p_a.cutoff_R);
  p_a.psi = tbl["psi"].value_or(p_a.psi);
}

} // namespace eonc::pot
