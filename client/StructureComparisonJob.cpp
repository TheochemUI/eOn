
#include "StructureComparisonJob.h"
#include "Optimizer.h"
#include "Log.h"
#include "Matter.h"
#include "HelperFunctions.h"

StructureComparisonJob::StructureComparisonJob(Parameters *params)
{
    parameters = params;
}

StructureComparisonJob::~StructureComparisonJob(){ }

std::vector<std::string> StructureComparisonJob::run(void)
{
    std::vector<std::string> returnFiles;
    
    Matter *matter1 = new Matter(parameters);
    matter1->con2matter("matter1.con");

    return returnFiles;
}
