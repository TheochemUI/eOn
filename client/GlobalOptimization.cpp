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
#include "eon/GlobalOptimization.h"
#include <cstdio>

#include "eon/EonLogger.h"

GlobalOptimization::GlobalOptimization(const Parameters &params)
    : parameters{params} {}

GlobalOptimization::~GlobalOptimization() = default;

void GlobalOptimization::run() { EONC_LOG_INFO("HELLO FROM GO\n"); }
