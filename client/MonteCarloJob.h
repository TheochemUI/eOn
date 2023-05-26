#ifndef MONTECARLOJOB_H
#define MONTECARLOJOB_H

#include "Job.h"
#include "Parameters.h"

class MonteCarloJob : public Job {
public:
  MonteCarloJob(std::unique_ptr<Parameters> parameters)
      : Job(std::move(parameters)) {}
  ~MonteCarloJob(void) = default;
  std::vector<std::string> run(void);
};

#endif
