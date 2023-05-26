#include "GlobalOptimization.h"
#include "Log.h"
#include <stdio.h>

GlobalOptimization::GlobalOptimization(Parameters *params) {
  parameters = params;
}

GlobalOptimization::~GlobalOptimization(void) {}

void GlobalOptimization::run(void) { log("HELLO FROM GO\n"); }
