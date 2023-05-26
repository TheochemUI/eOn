
#include "StructureComparisonJob.h"
#include "Optimizer.h"
#include "Log.h"
#include "Matter.h"
#include "HelperFunctions.h"

std::vector<std::string> StructureComparisonJob::run(void)
{
    std::vector<std::string> returnFiles;
    
    Matter *matter1 = new Matter(params);
    matter1->con2matter("matter1.con");

    return returnFiles;
}
