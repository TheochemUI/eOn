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

#include "C_Structs.h"

namespace eonc::pot {
// Typically this is done by the caller, however it is here as a sanity check
void zeroForceOut(const size_t &nAtoms, ForceOut *efvd);
} // namespace eonc::pot
