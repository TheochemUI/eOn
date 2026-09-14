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
  }
  ~VASP() { cleanMemory(); }
  void initialize() {};
  void cleanMemory(void);
  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box);
  /// force() writes POSCAR, signals through NEWCAR and waits on FU, all at
  /// fixed names in the working directory, so two threads in it swap results.
  [[nodiscard]] bool isThreadSafe() const noexcept override { return false; }
  /// Every instance drives the one VASP process named by the static vaspPID
  /// through the one set of files, so a second instance buys no parallelism.
  [[nodiscard]] bool needsPerImageInstance() const noexcept override {
    return false;
  }
  //!< Delete the results and restart files VASP writes into the working
  //!< directory, so a calculation starts from an empty one.
  //!<
  //!< Destructive, and nothing in the client calls it: a caller invokes it
  //!< when it means to discard an earlier calculation's output. A caller
  //!< that means to restart from WAVECAR or CHGCAR must not, and neither
  //!< must one that wants the earlier OUTCAR. The
  //!< examples/akmc-vasp-slurm scripts do the same removals in shell, at the
  //!< point where the user asks for them.
  static void removeStaleFiles();

private:
  void writePOSCAR(long N, const double *R, const int *atomicNrs,
                   const double *box);
  void readFU(long N, double *F, double *U);
  void spawnVASP();
  bool vaspRunning();
  //!< Remove the files eOn and VASP signal each other through, so a run does
  //!< not inherit the handshake of one that ended in this directory.
  static void clearHandshakeFiles();
  static bool firstRun;
  static long vaspRunCount;
  static pid_t vaspPID;
};
