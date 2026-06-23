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
// serves as an interface between LAMMPS potentials maintained by SANDIA

#pragma once

#include "../../Parameters.h"
#include "../../Potential.h"

class LAMMPSPot : public Potential {

public:
  [[nodiscard]] bool needsPerImageInstance() const noexcept override {
    return true;
  }
  LAMMPSPot(const Parameters &p);
  ~LAMMPSPot();
  void cleanMemory();
  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box) override;

private:
  int lammpsThr{0};
#ifdef EONMPI
  MPI_Comm mpiComm;
#endif
  long numberOfAtoms{0};
  double oldBox[9]{};
  void *LAMMPSObj{nullptr};
  void makeNewLAMMPS(long N, const double *R, const int *atomicNrs,
                     const double *box);
  bool realunits{false};

#if !defined(EONMPI) && !defined(IS_WINDOWS)
  // Process-per-image evaluation.  NEB drives intermediate images on separate
  // std::threads; if each thread opened LAMMPS in this process they would all
  // share one MPI_COMM_WORLD and their concurrent reduction collectives would
  // collide (heap corruption / MPI_ERR_OP).  Instead each LAMMPSPot forks a
  // dedicated worker process that owns its LAMMPS instance, so every image runs
  // in its own process with its own MPI_COMM_WORLD and true parallelism.
  // Not available on Windows (no fork/pipe).
  int workerPid{-1};
  int reqFd{-1}; // parent writes requests here (child stdin side)
  int resFd{-1}; // parent reads results here (child stdout side)
  bool workerSpawned{false};

  // Fork the worker child on first use; child enters runWorkerLoop().
  void ensureWorker();
  // Child main loop: read requests, evaluate, write results; never returns.
  [[noreturn]] void runWorkerLoop();
  void stopWorker();
#endif
  // In-process LAMMPS force evaluation (used directly on Windows/MPI, and
  // inside the worker child on POSIX).
  void forceLocal(long N, const double *R, const int *atomicNrs, double *F,
                  double *U, const double *box);
};
