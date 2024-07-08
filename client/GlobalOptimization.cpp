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
#include "GlobalOptimization.h"
#include <stdio.h>

GlobalOptimization::GlobalOptimization(Parameters *params) {
  parameters = params;
}

GlobalOptimization::~GlobalOptimization(void) {}

void GlobalOptimization::run(void) { SPDLOG_INFO("HELLO FROM GO\n"); }
