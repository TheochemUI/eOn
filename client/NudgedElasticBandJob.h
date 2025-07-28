#ifndef NEBJOB_H
#define NEBJOB_H

#include "Job.h"
#include "Matter.h"
#include "NudgedElasticBand.h"
#include "Parameters.h"

class NudgedElasticBandJob : public Job {

public:
  NudgedElasticBandJob(std::unique_ptr<Parameters> parameters)
      : Job(std::move(parameters)), fCallsNEB{0} {
    m_log = spdlog::get("combi");
  }
  ~NudgedElasticBandJob(void) = default;
  std::vector<std::string> run(void);

private:
  // functions
  void printEndState(NudgedElasticBand::NEBStatus status);
  void saveData(NudgedElasticBand::NEBStatus status, NudgedElasticBand *neb);

  // variables
  std::vector<std::string> returnFiles;
  size_t fCallsNEB;
  std::shared_ptr<spdlog::logger> m_log;
};

#endif
