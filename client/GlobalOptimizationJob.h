#ifndef GLOBALOPTIMIZATION_HJOB_H
#define GLOBALOPTIMIZATION_HJOB_H

#include "Job.h"
#include "Matter.h"
#include "Parameters.h"

class GlobalOptimizationJob : public Job {
public:
  GlobalOptimizationJob(std::unique_ptr<Parameters> parameters)
      : Job(std::move(parameters)), nlmin{0}, ediff{1.E-1}, ekin{5.E-2},
        beta1{params->globalOptimizationBeta},
        beta2{params->globalOptimizationBeta},
        beta3{1. / params->globalOptimizationBeta},
        alpha1{1. / params->globalOptimizationAlpha},
        alpha2{params->globalOptimizationAlpha},
        mdmin{params->globalOptimizationMdmin}, fcallsMove{0}, firstStep{true},
        fcallsRelax{0}, monfile{fopen("monitoring.dat", "w")}, earrfile{fopen(
                                                                   "earr.dat",
                                                                   "w")} {

    log = spdlog::get("combi");
  }
  // etoler = parameters->globalOptimizationEtoler;
  // decisionMethod = "NPEW";
  ~GlobalOptimizationJob(void) {
    fclose(monfile);
    fclose(earrfile);
  };
  void hoppingStep(long, Matter *, Matter *);
  void decisionStep(Matter *, Matter *);
  void report(Matter *);
  void acceptRejectNPEW(Matter *, Matter *);
  void acceptRejectBoltzmann(Matter *, Matter *);
  void analyze(Matter *, Matter *);
  void examineEscape(Matter *, Matter *);
  void applyMoveFeedbackMD(void);
  void applyDecisionFeedback(void);
  void mdescape(Matter *);
  void randomMove(Matter *);
  void insert(Matter *);
  size_t hunt(double);
  void velopt(Matter *);
  std::vector<std::string> run(void);
  double beta1;
  double beta2;
  double beta3;
  double alpha1;
  double alpha2;
  long mdmin;

private:
  size_t nlmin;
  long fcallsMove;
  long fcallsRelax;
  double ediff;
  double ekin;
  bool firstStep;
  std::vector<double> earr;
  string escapeResult;
  // string moveFeedbackMethod;
  // string decisionMethod;
  string decisionResult;
  string hoppingResult;
  FILE *monfile;
  FILE *earrfile;
  std::shared_ptr<spdlog::logger> log;
};

#endif
