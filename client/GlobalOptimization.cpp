/*
** This file is part of eON.
**
** SPDX-License-Identifier: BSD-3-Clause
**
** Copyright (c) 2010--present, eON Development Team
** All rights reserved.
**
** Repo:
** https://github.com/TheochemUI/eON
*/
#include "GlobalOptimization.h"

GlobalOptimization::GlobalOptimization(Parameters *params) {
  parameters = params;
}

GlobalOptimization::~GlobalOptimization(void) {}

void GlobalOptimization::run(void) { SPDLOG_INFO("HELLO FROM GO\n"); }
