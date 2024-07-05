#pragma once

#include "Job.h"
#include "Parameters.h"

class PrefactorJob : public Job {
public:
  PrefactorJob(std::unique_ptr<Parameters> parameters)
      : Job(std::move(parameters)) {}
  ~PrefactorJob(void) = default;
  std::vector<std::string> run(void);
  // Ugly but OK for now I guess
  static const char PREFACTOR_REACTANT[];
  static const char PREFACTOR_SADDLE[];
  static const char PREFACTOR_PRODUCT[];
};
