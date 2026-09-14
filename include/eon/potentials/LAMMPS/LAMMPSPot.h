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

#include "eon/Parameters.h"
#include "eon/Potential.h"

#include <mutex>

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
  // Respawns allowed after a worker times out, dies, or reports an error.
  // A single transient failure must not poison the job: every later
  // evaluation would return the impassable wall, no minimisation could ever
  // meet its force criterion, and the search would be discarded as a minimum
  // that failed to converge. Respawning is safe now that the teardown is
  // bounded rather than waiting on a wedged child forever. The budget is
  // finite so a worker that cannot be revived still ends the search instead
  // of looping.
  int workerRespawnsLeft{3};
  // Serialises the request/response exchange with the worker. eOn minimises
  // the two endpoints of a saddle concurrently, and when both share this
  // instance the two threads interleave writes and reads on the same pipe.
  // The protocol is a bare byte stream with no framing, so an interleaved
  // exchange is read as corrupt: the worker reports an evaluation error, the
  // next send finds a closed pipe, and the worker dies, all within the first
  // three force calls of the minimisation. Uncontended when instances really
  // are per-image.
  std::mutex workerMutex;
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
