/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
**
** Copyright (c) 2010--present, eOn Development Team
** All rights reserved.
**
** Repo:
** https://github.com/TheochemUI/eOn
*/
#pragma once

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
