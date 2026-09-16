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

/// Library surface for embedding eOn (IDPP/SIDPP paths, NEB, dimer).
/// Implementation headers stay under include/eon/; include this from
/// an external project instead of reaching into Job/*.h.
///
/// ABI (eOn-86bs, narrow cut):
/// - `Parameters` ships a private `Impl` for load state (`last_load_source`,
///   `last_load_error`). Option-group structs stay in the installed header,
///   so `sizeof(Parameters)` is still not ABI-stable.
/// - `Parameters.h` does not include Eigen. `Matter` still does: public
///   accessors return Eigen types, and hiding those members is a later cut.

#include "HelperFunctions.h"
#include "ImprovedDimer.h"
#include "JobResult.h"
#include "Lanczos.h"
#include "Matter.h"
#include "MinModeSaddleSearch.h"
#include "NEBInitialPaths.hpp"
#include "NudgedElasticBand.h"
#include "Optimizer.h"
#include "Parameters.h"
#include "PotCapabilities.h"
#include "Potential.h"
