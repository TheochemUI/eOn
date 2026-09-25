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
/// ABI:
/// - `Parameters` stores load state and option groups in a private `Impl`.
///   `sizeof(Parameters)` is that pointer. Const accessors and
///   `ParametersLoadAccess` are the read and write surface. Option-group
///   types remain in the installed `ParametersOptions.h`.
/// - `Parameters.h` does not include Eigen. `Matter` still does, because
///   accessors return Eigen types. Those accessors are defined in the
///   library. The matrices live in `Matter::Impl`, so `sizeof(Matter)` does
///   not embed them.

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
