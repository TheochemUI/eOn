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

#include "eon/Potential.h"

class VASP : public Potential {

public:
  VASP(const Parameters &p)
      : Potential(p) {
    vaspRunCount++;
    removeStaleFiles();
  }
  ~VASP() { cleanMemory(); }
  void initialize() {};
  void cleanMemory(void);
  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box);
  //!< Delete VASP output left in the working directory by an earlier run.
  //!< Called from the constructor, so constructing a VASP potential in a
  //!< directory that holds results deletes them.
  static void removeStaleFiles();

private:
  void writePOSCAR(long N, const double *R, const int *atomicNrs,
                   const double *box);
  void readFU(long N, double *F, double *U);
  void spawnVASP();
  bool vaspRunning();
  static bool firstRun;
  static long vaspRunCount;
  static pid_t vaspPID;
};
