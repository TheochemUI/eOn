#ifndef SAFEHYPERJOB_H
#define SAFEHYPERJOB_H

#include "Job.h"
#include "Matter.h"
#include "Parameters.h"

class SafeHyperJob : public Job {
public:
  SafeHyperJob(std::unique_ptr<Parameters> parameters)
      : Job(std::move(parameters)) {
    log = spdlog::get("combi");
  }
  ~SafeHyperJob(void) = default;
  std::vector<std::string> run(void);

private:
  shared_ptr<spdlog::logger> log;
  int dynamics();
  long refine(Matter *mdBuffer[], long length, Matter *reactant);
  bool checkState(Matter *current, Matter *reactant);
  void saveData(int status);
  void dephase();

  Matter *current;
  Matter *reactant;
  Matter *saddle;
  Matter *final_img;
  Matter *final_img_tmp;
  Matter *product;

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
  double transitionPot;
  double *timeBuffer;
  double *biasBuffer;

  std::vector<std::string> returnFiles;
};

#endif
