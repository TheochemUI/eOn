#pragma once

#include "client/potentials/LJ/LJ.h"
#include "client/potentials/Morse/Morse.h"

#ifdef WITH_ASE_ORCA
#include "client/potentials/ASE_ORCA/ASE_ORCA.h"
#endif

namespace eonc::pot {
void from_toml(LJ::Params &, const toml::node_view<const toml::node> &);
void from_toml(Morse::Params &, const toml::node_view<const toml::node> &);

#ifdef WITH_ASE_ORCA
void from_toml(ASEOrcaPot::Params &, const toml::node_view<const toml::node> &);
#endif
} // namespace eonc::pot
