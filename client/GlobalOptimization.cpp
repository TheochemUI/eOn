#include "GlobalOptimization.h"

GlobalOptimization::GlobalOptimization(Parameters *params) {
  parameters = params;
}

GlobalOptimization::~GlobalOptimization(void) {}

void GlobalOptimization::run(void) { SPDLOG_INFO("HELLO FROM GO\n"); }
