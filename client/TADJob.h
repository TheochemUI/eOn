
#ifndef TADJOB_H
#define TADJOB_H

#include "Job.h"
#include "MinModeSaddleSearch.h"
#include "Parameters.h"

class TADJob : public Job {
public:
  TADJob(std::unique_ptr<Parameters> parameters) : Job(std::move(parameters)) {
    log = spdlog::get("combi");
  }
  ~TADJob() = default;
  std::vector<std::string> run(void);

private:
  int dynamics();
  long refine(std::vector<std::shared_ptr<Matter>> buff, long length,
              Matter *reactant);
  bool checkState(Matter *current, Matter *reactant);
  void saveData(int status);
  void dephase();
  bool saddleSearch(std::shared_ptr<Matter> cross);

  std::shared_ptr<Matter> current;
  std::shared_ptr<Matter> reactant;
  std::shared_ptr<Matter> saddle;
  std::shared_ptr<Matter> crossing;
  std::shared_ptr<Matter> final_state;
  std::shared_ptr<Matter> final_tmp;
  std::shared_ptr<Matter> product;
  std::shared_ptr<spdlog::logger> log;

  bool metaStateFlag;
  bool newStateFlag;

  long minimizeFCalls;
  long mdFCalls;
  long dephaseFCalls;
  long refineFCalls;

  long transitionStep;

  double time;
  double minCorrectedTime;
  double transitionTime;
  double barrier;
  double *timeBuffer;
  MinModeSaddleSearch *dimerSearch;

  std::vector<std::string> returnFiles;
};

#endif
