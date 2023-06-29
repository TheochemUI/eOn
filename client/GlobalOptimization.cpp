#include "GlobalOptimization.h"
#include <stdio.h>

GlobalOptimization::GlobalOptimization(Parameters *params) {
  parameters = params;
}

GlobalOptimization::~GlobalOptimization(void) {}

void GlobalOptimization::run(void) { SPDLOG_INFO("HELLO FROM GO\n"); }
