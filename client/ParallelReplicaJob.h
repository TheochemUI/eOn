#ifndef PARALLELREPLICAJOB_H
#define PARALLELREPLICAJOB_H

#include "Job.h"
#include "Parameters.h"

class ParallelReplicaJob : public Job {
public:
  ParallelReplicaJob(std::unique_ptr<Parameters> parameters)
      : Job(std::move(parameters)) {
    log = spdlog::get("combi");
  }
  ~ParallelReplicaJob(void) = default;
  std::vector<std::string> run(void);

private:
  std::vector<std::string> returnFiles;
  Matter *reactant;

  void dephase(Matter *trajectory);
  int refineTransition(std::vector<Matter *> MDSnapshots, bool fake = false);
  std::shared_ptr<spdlog::logger> log;
};

#endif
