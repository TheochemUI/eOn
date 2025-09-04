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

#include "../../Potential.h"

class VASP : public Potential {

public:
  VASP(shared_ptr<Parameters> p)
      : Potential(p) {
    vaspRunCount++;
    // deleting leftovers from previous run
    system("rm -f TMPCAR");
    system("rm -f CHG");
    system("rm -f CHGCAR");
    system("rm -f CONTCAR");
    system("rm -f DOSCAR");
    system("rm -f EIGENVAL");
    system("rm -f IBZKPT");
    system("rm -f NEWCAR");
    system("rm -f FU");
    system("rm -f OSZICAR");
    system("rm -f OUTCAR");
    system("rm -f PCDAT");
    system("rm -f POSCAR");
    system("rm -f TMPCAR");
    system("rm -f WAVECAR");
    system("rm -f XDATCAR");
  }
  ~VASP() { cleanMemory(); }
  void initialize() {};
  void cleanMemory(void);
  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box);

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
