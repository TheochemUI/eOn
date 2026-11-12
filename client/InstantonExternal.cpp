/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
**
** Copyright (c) 2010--present, eOn Development Team
** All rights reserved.
*/
#include "eon/InstantonExternal.h"

namespace eonc {
namespace {
ExternalInstanton g_external{};
}

void set_external_instanton(ExternalInstanton hook) { g_external = hook; }

ExternalInstanton external_instanton() { return g_external; }

} // namespace eonc
