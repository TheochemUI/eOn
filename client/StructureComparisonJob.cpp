#include "StructureComparisonJob.h"
#include "HelperFunctions.h"
#include "Matter.h"
#include "Optimizer.h"

std::vector<std::string> StructureComparisonJob::run(void) {
  std::vector<std::string> returnFiles;

  auto matter1 = std::make_unique<Matter>(pot, params);
  matter1->con2matter("matter1.con");

  return returnFiles;
}
