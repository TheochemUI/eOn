
#ifndef TESTJOB_H
#define TESTJOB_H

#include "ConjugateGradients.h"
#include "Job.h"
#include "Parameters.h"

class TestJob : public Job {
public:
  TestJob(Parameters *params);
  ~TestJob(void);
  std::vector<std::string> run(void);

private:
  Parameters *parameters;
  double tolerance;
  void checkFullSearch(void);
  void checkMode(void);
  void checkPotentials(void);
  double getEnergyDiff(string potTag, double refEnergy);
  double getForceDiff(string potTag, double refForce);
};

#endif
