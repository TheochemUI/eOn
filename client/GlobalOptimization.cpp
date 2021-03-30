#include "GlobalOptimization.h"
#include <stdio.h>
#include "Log.h"

GlobalOptimization::GlobalOptimization(Parameters *params)
{
    parameters = params;
}

GlobalOptimization::~GlobalOptimization(void)
{

}

void GlobalOptimization::run(void)
{
    log("HELLO FROM GO\n");
}
