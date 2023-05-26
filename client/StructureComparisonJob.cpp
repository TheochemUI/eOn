
#include "StructureComparisonJob.h"
#include "HelperFunctions.h"
#include "Log.h"
#include "Matter.h"
#include "Optimizer.h"

std::vector<std::string> StructureComparisonJob::run(void) {
  std::vector<std::string> returnFiles;

  Matter *matter1 = new Matter(params);
  matter1->con2matter("matter1.con");

  return returnFiles;
}
