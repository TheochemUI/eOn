#pragma once
#include "Job.h"
#include "Matter.h"
#include "Parameters.h"

class GlobalOptimizationJob : public Job {
public:
  GlobalOptimizationJob(std::unique_ptr<Parameters> parameters)
      : Job(std::move(parameters)),
        beta1{params->globopt.beta},
        beta2{params->globopt.beta},
        beta3{1. / params->globopt.beta},
        alpha1{1. / params->globopt.alpha},
        alpha2{params->globopt.alpha},
        mdmin{params->globopt.mdmin},
        nlmin{0},
        fcallsMove{0},
        fcallsRelax{0},
        ediff{1.E-1},
        ekin{5.E-2},
        firstStep{true},
        monfile{fopen("monitoring.dat", "w")},
        earrfile{fopen("earr.dat", "w")} {

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
  std::string escapeResult;
  // string moveFeedbackMethod;
  // string decisionMethod;
  std::string decisionResult;
  std::string hoppingResult;
  FILE *monfile;
  FILE *earrfile;
  std::shared_ptr<spdlog::logger> log;
};
