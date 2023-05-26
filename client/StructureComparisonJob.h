
#ifndef STRUCTURECOMPARISONJOB_H
#define STRUCTURECOMPARISONJOB_H

#include "Job.h"
#include "Parameters.h"

class StructureComparisonJob : public Job {
public:
  StructureComparisonJob(std::unique_ptr<Parameters> parameters)
      : Job(std::move(parameters)) {}
  ~StructureComparisonJob(void) = default;
  std::vector<std::string> run(void);
};

#endif
