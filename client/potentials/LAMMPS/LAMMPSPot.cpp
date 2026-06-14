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
#include "LAMMPSPot.h"
#include "LammpsLoader.h"

#include <cstring>
#include <filesystem>
#include <format>
#include <fstream>
#include <map>
#include <string>

#ifndef EONMPI
#include <cerrno>
#include <cstdlib>
#include <vector>

#include <sys/wait.h>
#include <unistd.h>
#endif

#ifdef EONMPI
#define LAMMPS_LIB_MPI
#endif

LAMMPSPot::LAMMPSPot(const Parameters &p)
    : Potential(p),
      lammpsThr{p.potential_options.LAMMPSThreads}
#ifdef EONMPI
      ,
      mpiComm{p.potential_options.MPIClientComm}
#endif
{
  // Fail fast if LAMMPS library not available
  eonc::LammpsLoader::instance().require_loaded();
#ifndef EONMPI
  // Fork the worker NOW, at construction, before this process ever opens a
  // LAMMPS instance (and thus before liblammps initialises MPI).  Open MPI
  // does not support using MPI in a process that called MPI_Init before fork,
  // so the worker must be spawned from a still-MPI-clean parent.  Every
  // LAMMPSPot -- endpoints and per-image alike -- runs its LAMMPS in its own
  // child process, so the parent never initialises MPI at all.
  ensureWorker();
#endif
}

LAMMPSPot::~LAMMPSPot() { cleanMemory(); }

void LAMMPSPot::cleanMemory() {
#ifndef EONMPI
  stopWorker();
#endif
  if (LAMMPSObj != nullptr) {
    eonc::LammpsLoader::instance().close(LAMMPSObj);
    LAMMPSObj = nullptr;
  }
}

#ifndef EONMPI
// ---------------------------------------------------------------------------
// Process-per-image worker plumbing
// ---------------------------------------------------------------------------
namespace {
// Blocking read/write of exactly n bytes over a pipe.  Returns false on EOF or
// error, so a dead peer is detected rather than silently producing garbage.
bool readExact(int fd, void *buf, size_t n) {
  auto *p = static_cast<char *>(buf);
  while (n > 0) {
    ssize_t r = read(fd, p, n);
    if (r <= 0) {
      if (r < 0 && errno == EINTR)
        continue;
      return false;
    }
    p += r;
    n -= static_cast<size_t>(r);
  }
  return true;
}
bool writeExact(int fd, const void *buf, size_t n) {
  const auto *p = static_cast<const char *>(buf);
  while (n > 0) {
    ssize_t w = write(fd, p, n);
    if (w < 0) {
      if (errno == EINTR)
        continue;
      return false;
    }
    p += w;
    n -= static_cast<size_t>(w);
  }
  return true;
}
} // namespace

void LAMMPSPot::ensureWorker() {
  if (workerSpawned)
    return;

  int reqPipe[2]; // parent -> child
  int resPipe[2]; // child -> parent
  if (pipe(reqPipe) != 0 || pipe(resPipe) != 0) {
    throw std::runtime_error("LAMMPSPot: failed to create worker pipes");
  }

  // Fork BEFORE opening any LAMMPS instance in this process, so MPI is first
  // initialised inside the child.  Each child is its own process with its own
  // MPI_COMM_WORLD; concurrent children never share a communicator.
  pid_t pid = fork();
  if (pid < 0) {
    throw std::runtime_error("LAMMPSPot: fork for worker failed");
  }

  if (pid == 0) {
    // Child: keep reqPipe read end and resPipe write end.
    close(reqPipe[1]);
    close(resPipe[0]);
    reqFd = reqPipe[0];
    resFd = resPipe[1];
    runWorkerLoop(); // never returns
  }

  // Parent: keep reqPipe write end and resPipe read end.
  close(reqPipe[0]);
  close(resPipe[1]);
  reqFd = reqPipe[1];
  resFd = resPipe[0];
  workerPid = pid;
  workerSpawned = true;
}

void LAMMPSPot::runWorkerLoop() {
  // Running in the forked child.  Evaluate forces with an in-process LAMMPS
  // (this child's own MPI_COMM_WORLD) and stream results back to the parent.
  for (;;) {
    long N = 0;
    if (!readExact(reqFd, &N, sizeof(N))) {
      _exit(0); // request pipe closed -> shut down cleanly
    }
    if (N < 0) {
      _exit(0); // explicit shutdown sentinel from stopWorker()
    }
    std::vector<int> atomicNrs(static_cast<size_t>(N));
    std::vector<double> R(static_cast<size_t>(3 * N));
    double box[9];
    if (!readExact(reqFd, atomicNrs.data(), sizeof(int) * static_cast<size_t>(N)) ||
        !readExact(reqFd, box, sizeof(box)) ||
        !readExact(reqFd, R.data(), sizeof(double) * static_cast<size_t>(3 * N))) {
      _exit(1);
    }

    std::vector<double> F(static_cast<size_t>(3 * N), 0.0);
    double U = 0.0;
    int status = 0;
    try {
      forceLocal(N, R.data(), atomicNrs.data(), F.data(), &U, box);
    } catch (...) {
      status = 1;
    }

    if (!writeExact(resFd, &status, sizeof(status)) ||
        !writeExact(resFd, &U, sizeof(U)) ||
        !writeExact(resFd, F.data(), sizeof(double) * static_cast<size_t>(3 * N))) {
      _exit(1);
    }
  }
}

void LAMMPSPot::stopWorker() {
  if (!workerSpawned)
    return;
  if (reqFd >= 0) {
    // Send an explicit shutdown sentinel, then close.  A sentinel (rather than
    // relying on pipe EOF) guarantees the child exits even when sibling worker
    // processes hold an inherited copy of this write end.
    long sentinel = -1;
    writeExact(reqFd, &sentinel, sizeof(sentinel));
    close(reqFd);
    reqFd = -1;
  }
  if (resFd >= 0) {
    close(resFd);
    resFd = -1;
  }
  if (workerPid > 0) {
    int st = 0;
    waitpid(workerPid, &st, 0);
    workerPid = -1;
  }
  workerSpawned = false;
}
#endif // !EONMPI

void LAMMPSPot::force(long N, const double *R, const int *atomicNrs, double *F,
                      double *U, double *variance, const double *box) {
  variance = nullptr;

#ifdef EONMPI
  forceLocal(N, R, atomicNrs, F, U, box);
#else
  // Drive the dedicated worker process so this image's LAMMPS runs in its own
  // process (own MPI_COMM_WORLD).  Per-image NEB threads thus evaluate forces
  // as truly concurrent processes with no shared-communicator contention.
  ensureWorker();

  if (!writeExact(reqFd, &N, sizeof(N)) ||
      !writeExact(reqFd, atomicNrs, sizeof(int) * static_cast<size_t>(N)) ||
      !writeExact(reqFd, box, sizeof(double) * 9) ||
      !writeExact(reqFd, R, sizeof(double) * static_cast<size_t>(3 * N))) {
    throw std::runtime_error("LAMMPSPot: failed to send request to worker");
  }

  int status = 0;
  if (!readExact(resFd, &status, sizeof(status)) ||
      !readExact(resFd, U, sizeof(double)) ||
      !readExact(resFd, F, sizeof(double) * static_cast<size_t>(3 * N))) {
    throw std::runtime_error("LAMMPSPot: worker process died during force eval");
  }
  if (status != 0) {
    throw std::runtime_error("LAMMPSPot: worker reported a force evaluation error");
  }
#endif
}

void LAMMPSPot::forceLocal(long N, const double *R, const int *atomicNrs,
                           double *F, double *U, const double *box) {
  auto &lmp = eonc::LammpsLoader::instance();

  bool newLammps = false;
  for (int i = 0; i < 9; i++) {
    if (oldBox[i] != box[i])
      newLammps = true;
  }
  if (numberOfAtoms != N)
    newLammps = true;
  if (newLammps) {
    makeNewLAMMPS(N, R, atomicNrs, box);
  }
  if (!LAMMPSObj) {
    throw std::runtime_error("Should have a LAMMPS instance by now");
  }

  lmp.scatter_atoms(LAMMPSObj, "x", 1, 3, const_cast<double *>(R));
  lmp.command(LAMMPSObj, "run 1 pre no post no");

  auto *pe =
      static_cast<double *>(lmp.extract_variable(LAMMPSObj, "pe", nullptr));
  *U = *pe;
  free(pe);

  auto *fx =
      static_cast<double *>(lmp.extract_variable(LAMMPSObj, "fx", "all"));
  auto *fy =
      static_cast<double *>(lmp.extract_variable(LAMMPSObj, "fy", "all"));
  auto *fz =
      static_cast<double *>(lmp.extract_variable(LAMMPSObj, "fz", "all"));

  for (long i = 0; i < N; i++) {
    F[3 * i + 0] = fx[i];
    F[3 * i + 1] = fy[i];
    F[3 * i + 2] = fz[i];
  }

  // Convert kCal/mol -> eV if LAMMPS is using real units
  if (realunits) {
    constexpr double kcalPerEv = 23.0609;
    *U /= kcalPerEv;
    for (long i = 0; i < 3 * N; i++) {
      F[i] /= kcalPerEv;
    }
  }

  free(fx);
  free(fy);
  free(fz);
}

void LAMMPSPot::makeNewLAMMPS(long N, const double *R, const int *atomicNrs,
                              const double *box) {
  auto &lmp = eonc::LammpsLoader::instance();

  numberOfAtoms = N;
  std::memcpy(oldBox, box, 9 * sizeof(double));

  if (LAMMPSObj != nullptr) {
    eonc::LammpsLoader::instance().close(LAMMPSObj);
    LAMMPSObj = nullptr;
  }

  // Map atomic numbers to LAMMPS type indices (1-based)
  std::map<int, int> type_map;
  int ntypes = 0;
  for (long i = 0; i < N; i++) {
    if (type_map.count(atomicNrs[i]) == 0) {
      type_map.insert({atomicNrs[i], ++ntypes});
    }
  }

#ifdef EONMPI
  const char *lmpargv[] = {"liblammps", "-log", "none",    "-echo", "log",
                           "-screen",   "none", "-suffix", "omp"};
  int lmpargc = sizeof(lmpargv) / sizeof(const char *);
  if (!lmp.open_mpi) {
    throw std::runtime_error(
        "LAMMPS library found but lacks MPI support (lammps_open not found).\n"
        "Install an MPI-enabled LAMMPS build.");
  }
  MPI_Comm inst_comm = MPI_COMM_NULL;
  MPI_Comm_dup(mpiComm, &inst_comm);  // private comm per per-image instance
  LAMMPSObj =
      lmp.open_mpi(lmpargc, const_cast<char **>(lmpargv), inst_comm, nullptr);
#else
  const char *lmpargv[] = {"liblammps", "-log",    "none", "-echo",
                           "log",       "-screen", "none"};
  int lmpargc = sizeof(lmpargv) / sizeof(const char *);
  LAMMPSObj = lmp.open_no_mpi(lmpargc, const_cast<char **>(lmpargv), nullptr);
#endif

  if (lammpsThr > 0) {
    std::string cmd = std::format("package omp {} force/neigh", lammpsThr);
    lmp.command(LAMMPSObj, cmd.c_str());
  }

  // Detect units from in.lammps: look for "#!units real" marker
  realunits = false;
  if (std::filesystem::exists("in.lammps")) {
    std::ifstream infile("in.lammps");
    std::string line;
    while (std::getline(infile, line)) {
      if (line == "#!units real") {
        realunits = true;
        break;
      }
    }
  } else {
    EONC_LOG_ERROR("[LAMMPS] in.lammps not found in working directory");
    return;
  }

  if (realunits) {
    lmp.command(LAMMPSObj, "units real");
  } else {
    lmp.command(LAMMPSObj, "units metal");
  }

  lmp.command(LAMMPSObj, "atom_style charge");
  lmp.command(LAMMPSObj, "atom_modify map array sort 0 0");
  lmp.command(LAMMPSObj, "neigh_modify delay 1");

  // Define periodic cell (prism for non-orthorhombic)
  std::string region_cmd =
      std::format("region cell prism 0 {} 0 {} 0 {} {} {} {} units box", box[0],
                  box[4], box[8], box[3], box[6], box[7]);
  lmp.command(LAMMPSObj, region_cmd.c_str());

  std::string create_box_cmd = std::format("create_box {} cell", ntypes);
  lmp.command(LAMMPSObj, create_box_cmd.c_str());

  // Initialize atoms
  for (long i = 0; i < N; i++) {
    std::string atom_cmd =
        std::format("create_atoms {} single {} {} {} units box",
                    type_map[atomicNrs[i]], 0.0, 0.0, 0.0);
    lmp.command(LAMMPSObj, atom_cmd.c_str());
  }

  lmp.command(LAMMPSObj, "mass * 1.0");

  // Load user LAMMPS input script
  lmp.file(LAMMPSObj, "in.lammps");

  // Define variables for force/energy extraction
  lmp.command(LAMMPSObj, "variable fx atom fx");
  lmp.command(LAMMPSObj, "variable fy atom fy");
  lmp.command(LAMMPSObj, "variable fz atom fz");
  lmp.command(LAMMPSObj, "variable pe equal pe");
}
