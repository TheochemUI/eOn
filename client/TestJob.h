#pragma once

#include "Job.h"
#include "Parameters.h"

class TestJob : public Job {
public:
  TestJob(std::unique_ptr<Parameters> params)
      : Job(std::move(params)),
        tolerance{0.01} {}
  ~TestJob(void) = default;
  std::vector<std::string> run(void);

private:
  double tolerance;
  void checkFullSearch(void);
  void checkMode(void);
  void checkPotentials(void);
  double getEnergyDiff(std::string potTag, double refEnergy);
  double getForceDiff(std::string potTag, double refForce);
};
