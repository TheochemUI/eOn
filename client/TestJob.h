
#ifndef TESTJOB_H
#define TESTJOB_H

#include "ConjugateGradients.h"
#include "Job.h"
#include "Parameters.h"

class TestJob : public Job {
public:
  TestJob(std::unique_ptr<Parameters> params)
      : Job(std::move(params)), tolerance{0.01} {}
  ~TestJob(void) = default;
  std::vector<std::string> run(void);

private:
  double tolerance;
  void checkFullSearch(void);
  void checkMode(void);
  void checkPotentials(void);
  double getEnergyDiff(string potTag, double refEnergy);
  double getForceDiff(string potTag, double refForce);
};

#endif
