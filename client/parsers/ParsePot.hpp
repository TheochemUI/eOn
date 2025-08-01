/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
**
** Copyright (c) 2010--present, eOn Development Team
** All rights reserved.
**
** Repo:
** https://github.com/TheochemUI/eOn
*/
#pragma once

#include "client/potentials/LJ/LJ.h"
#include "client/potentials/Morse/Morse.h"

#ifdef WITH_ASE_ORCA
#include "client/potentials/ASE_ORCA/ASE_ORCA.h"
#endif

#ifdef LAMMPS_POT
#include "client/potentials/LAMMPS/LAMMPSPot.h"
#endif

namespace eonc::pot {
void from_toml(LJ::Params &, const toml::node_view<const toml::node> &);
void from_toml(Morse::Params &, const toml::node_view<const toml::node> &);

#ifdef WITH_ASE_ORCA
void from_toml(ASEOrcaPot::Params &, const toml::node_view<const toml::node> &);
#endif

#ifdef LAMMPS_POT
void from_toml(LAMMPSPot::Params &, const toml::node_view<const toml::node> &);
#endif
} // namespace eonc::pot
