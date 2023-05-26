#ifndef HESSIANJOB_H
#define HESSIANJOB_H

#include "Job.h"
#include "Parameters.h"

class HessianJob : public Job {
public:
  HessianJob(std::unique_ptr<Parameters> parameters)
      : Job(std::move(parameters)) {}
  ~HessianJob(void) = default;
  std::vector<std::string> run(void);
};

#endif
