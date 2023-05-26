#ifndef MINIMAZATIONJOB_H
#define MINIMAZATIONJOB_H

#include "Job.h"
#include "Parameters.h"

class MinimizationJob : public Job {
public:
  MinimizationJob(std::unique_ptr<Parameters> parameters)
      : Job(std::move(parameters)), fcalls{0} {}
  ~MinimizationJob(void) = default;
  std::vector<std::string> run(void);

private:
  size_t fcalls;
  RunStatus status;
};

#endif
