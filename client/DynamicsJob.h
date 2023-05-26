#ifndef DYNAMICSJOB_H
#define DYNAMICSJOB_H

#include "Job.h"
#include "Matter.h"
#include "Parameters.h"

class DynamicsJob : public Job {

public:
  DynamicsJob(std::unique_ptr<Parameters> parameters)
      : Job(std::move(parameters)) {}
  ~DynamicsJob(void) = default;
  std::vector<std::string> run(void);
  std::vector<std::string> returnFiles;
};

#endif
