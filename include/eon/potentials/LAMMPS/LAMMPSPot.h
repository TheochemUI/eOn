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

#include <cerrno>
#include <cstdint>
#include <mutex>
#include <string>
#include <vector>

#ifdef EONMPI
#include <mpi.h>
#endif

namespace eonc {
/// True when waitpid has collected the child. EINTR is not a collection.
inline bool lammpsWorkerReaped(long got, long child, int err) {
  if (got == child) {
    return true;
  }
  return got < 0 && err != EINTR;
}

/// Byte offset at which to keep copying a LAMMPS screen file. A new LAMMPS
/// open truncates that file, so a restart or a shorter file reads from the
/// start.
inline std::int64_t lammpsScreenCursor(std::int64_t pos, std::int64_t fileSize,
                                       bool restarted) {
  if (restarted || pos < 0 || fileSize < pos) {
    return 0;
  }
  return pos;
}

/// LAMMPS argv. The log file stays off. A screen path is copied into the
/// process logger by the caller.
inline std::vector<std::string> lammpsOpenArgs(bool logging, bool with_omp,
                                               const std::string &screen) {
  std::vector<std::string> args{"liblammps", "-echo", "screen", "-log", "none"};
  if (logging) {
    args.insert(args.end(), {"-screen", screen});
  } else {
    args.insert(args.end(), {"-screen", "none"});
  }
  if (with_omp) {
    args.insert(args.end(), {"-suffix", "omp"});
  }
  return args;
}
} // namespace eonc

namespace eonc {
class ILammpsLoader;
}

class LAMMPSPot : public eonc::Potential {

public:
  [[nodiscard]] bool needsPerImageInstance() const noexcept override {
    return true;
  }
  /// Production: process-default LammpsLoader and POSIX worker isolation.
  explicit LAMMPSPot(const eonc::Parameters &p);
  /// Test seam: injected loader, no worker fork, no process-default load.
  LAMMPSPot(const eonc::Parameters &p, eonc::ILammpsLoader &loader);
  ~LAMMPSPot();
  void cleanMemory();
  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box) override;
  void setFixedMask(long nAtoms, const double *isFixed) override;

private:
  LAMMPSPot(const eonc::Parameters &p, eonc::ILammpsLoader &loader,
            bool isolate_worker);
  eonc::ILammpsLoader &loader_;
  int lammpsThr{0};
  bool lammpsLogging_{false};
  int lammpsLogIndex_{0};
  std::mutex maskMutex_;
  // client_lammps-N.log. The worker child writes it. The parent, after the
  // child has finished a force call, copies new lines into the process log.
  // The child does not touch that logger: quill does not survive fork.
  std::string lammpsScreenPath_;
  std::int64_t lammpsScreenPos_{0};
  // Consumed by drainLammpsScreen. A new LAMMPS open truncates the screen
  // file, so the next copy starts at the beginning.
  bool lammpsScreenRestart_{false};
  bool workerChild_{false};
#ifdef EONMPI
  MPI_Comm mpiComm;
#endif
  long numberOfAtoms{0};
  double oldBox[9]{};
  void *LAMMPSObj{nullptr};
  void makeNewLAMMPS(long N, const double *R, const int *atomicNrs,
                     const double *box);
  void applySetforce(long N);
  bool realunits{false};
  std::vector<double> fixedMask_;
  long maskN_{0};
  // Covers fixedMask_ and, on the in-process paths, LAMMPSObj. The worker
  // pipe uses the same mutex: a shared instance must not update the mask
  // while another thread copies it or evaluates a force.
  std::mutex workerMutex;

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
  int workerPid{-1};
  int reqFd{-1}; // parent writes requests here (child stdin side)
  int resFd{-1}; // parent reads results here (child stdout side)
  bool workerSpawned{false};
  // Geometry of the last request. The child rebuilds LAMMPS, and truncates
  // the screen file, when the atom count or the cell changes.
  bool screenHaveGeom_{false};
  long screenAtoms_{0};
  double screenBox_[9]{};

  // Fork the worker child on first use; child enters runWorkerLoop().
  void ensureWorker();
  bool noteScreenGeometry(long N, const double *box);
  // Child main loop: read requests, evaluate, write results; never returns.
  [[noreturn]] void runWorkerLoop();
  void stopWorker();
#endif
  // In-process LAMMPS force evaluation (used directly on Windows/MPI, and
  // inside the worker child on POSIX).
  void forceLocal(long N, const double *R, const int *atomicNrs, double *F,
                  double *U, const double *box);
  void lammpsCommand(const char *cmd);
  void drainLammpsScreen();
};
