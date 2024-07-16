#include "ParseTOML.hpp"

namespace eonc::pot {

void from_toml(LJ::Params &p_a, const toml::node_view<const toml::node> &tbl) {
  p_a.u0 = tbl["u0"].value_or(p_a.u0);
  p_a.cutoff_R = tbl["cutoff"].value_or(p_a.cutoff_R);
  p_a.psi = tbl["psi"].value_or(p_a.psi);
}

void from_toml(Morse::Params &p_a,
               const toml::node_view<const toml::node> &tbl) {
  p_a.De = tbl["Potential"]["Morse"]["De"].value_or(p_a.De);
  p_a.a = tbl["Potential"]["Morse"]["a"].value_or(p_a.a);
  p_a.re = tbl["Potential"]["Morse"]["re"].value_or(p_a.re);
  p_a.cutoff = tbl["Potential"]["Morse"]["cutoff"].value_or(p_a.cutoff);
}

#ifdef WITH_ASE_ORCA
void from_toml(ASEOrcaPot::Params &p_a,
               const toml::node_view<const toml::node> &tbl) {
  p_a.orca_path = tbl["ASE_ORCA"]["orca_path"].value_or(p_a.orca_path);
  p_a.orca_nproc = tbl["ASE_ORCA"]["nproc"].value_or(p_a.orca_nproc);
  p_a.simpleinput = tbl["ASE_ORCA"]["simpleinput"].value_or(p_a.simpleinput);
}
#endif
} // namespace eonc::pot
