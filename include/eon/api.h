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

#include "Matter.h"
#include "Parameters.h"
#include "Potential.h"
#include "PotCapabilities.h"
#include "JobResult.h"
#include "NEBInitialPaths.hpp"
#include "NudgedElasticBand.h"
#include "MinModeSaddleSearch.h"
#include "ImprovedDimer.h"
#include "Lanczos.h"
#include "Optimizer.h"
#include "HelperFunctions.h"
