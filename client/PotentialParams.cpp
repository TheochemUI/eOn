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
#include "eon/Parameters.h"
#include "eon/Potential.h"

namespace eonc {

Potential::Potential(PotType a_ptype, const Parameters &p)
    : Potential(a_ptype) {
  force_serial_ = !p.potential_options().thread_safe;
}

Potential::Potential(const Parameters &a_params)
    : Potential(a_params.potential_options().potential, a_params) {}

} // namespace eonc
