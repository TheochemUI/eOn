/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
**
** Copyright (c) 2010--present, eOn Development Team
** All rights reserved.
*/
#pragma once

namespace eonc {

/// External orbit. `evaluate` returns 0 and writes `action` when it owns
/// the rate. A null `evaluate` leaves the built-in ring polymer in place.
struct ExternalInstanton {
  int (*evaluate)(const double *position, int n_coord, const double *masses,
                  int n_atom, double temperature_k, double *action,
                  void *user);
  void *user;
};

void set_external_instanton(ExternalInstanton hook);
ExternalInstanton external_instanton();

} // namespace eonc
