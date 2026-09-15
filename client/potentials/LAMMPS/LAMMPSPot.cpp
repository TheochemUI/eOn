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
#include "eon/potentials/LAMMPS/LAMMPSPot.h"
#include "eon/EonLogger.h"
#include "eon/fpe_handler.h"
#include "eon/potentials/LAMMPS/LammpsLoader.h"

#include <cmath>
#include <cstring>
#include <filesystem>
#include <format>
#include <fstream>
#include <map>
#include <stdexcept>
#include <string>

#if !defined(EONMPI) && !defined(IS_WINDOWS)
#include <cerrno>
#include <cstdlib>
#include <vector>

#include <csignal>
#include <poll.h>
#include <sys/wait.h>
#include <unistd.h>
#endif

#ifdef EONMPI
#define LAMMPS_LIB_MPI
#endif

LAMMPSPot::LAMMPSPot(const eonc::Parameters &p)
    : eonc::Potential(p),
      lammpsThr{p.potential_options.LAMMPSThreads}
#ifdef EONMPI
      ,
      mpiComm{p.potential_options.MPIClientComm}
#endif
{
  // Fail fast if LAMMPS library not available
  eonc::LammpsLoader::instance().require_loaded();
#if !defined(EONMPI) && !defined(IS_WINDOWS)
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

void LAMMPSPot::setFixedMask(long nAtoms, const double *isFixed) {
  if (nAtoms <= 0 || isFixed == nullptr) {
    fixedMask_.clear();
    maskN_ = 0;
    return;
  }
  fixedMask_.assign(isFixed, isFixed + 3 * nAtoms);
  maskN_ = nAtoms;
}

void LAMMPSPot::applySetforce(long N) {
  if (LAMMPSObj == nullptr || maskN_ != N || fixedMask_.empty()) {
    return;
  }
  auto &lmp = eonc::LammpsLoader::instance();
  try {
    lmp.command(LAMMPSObj, "unfix eon_freeze");
  } catch (...) {
  }
  try {
    lmp.command(LAMMPSObj, "group eon_frozen delete");
  } catch (...) {
  }
  std::string ids;
  for (long i = 0; i < N; ++i) {
    if (fixedMask_[static_cast<size_t>(3 * i)] >= 0.5 &&
        fixedMask_[static_cast<size_t>(3 * i + 1)] >= 0.5 &&
        fixedMask_[static_cast<size_t>(3 * i + 2)] >= 0.5) {
      ids += std::format("{} ", i + 1);
    }
  }
  if (ids.empty()) {
    return;
  }
  lmp.command(LAMMPSObj, ("group eon_frozen id " + ids).c_str());
  lmp.command(LAMMPSObj, "fix eon_freeze eon_frozen setforce 0.0 0.0 0.0");
}

void LAMMPSPot::cleanMemory() {
#if !defined(EONMPI) && !defined(IS_WINDOWS)
  stopWorker();
#endif
  if (LAMMPSObj != nullptr) {
    eonc::LammpsLoader::instance().close(LAMMPSObj);
    LAMMPSObj = nullptr;
  }
}

#if !defined(EONMPI) && !defined(IS_WINDOWS)
// ---------------------------------------------------------------------------
// Process-per-image worker plumbing (POSIX only)
// ---------------------------------------------------------------------------
namespace {
// Report a geometry the worker could not evaluate as an impassable wall.
//
// Every failure path here used to throw, and nothing between the potential
// call and main catches it, so the client terminated. That loses the whole
// job, including searches that had already converged: a copper V/SIA search
// reached a saddle at 0.082 eV with a force of 0.027 eV/A and a curvature of
// -0.47, then died on the next evaluation before the result was written.
//
// A large finite energy with zeroed forces reads to the optimiser as a wall,
// so it backs out of the step and the search abandons this configuration and
// carries on. Callers stop the worker first; ensureWorker respawns it on the
// next evaluation.
void rejectGeometry(double *U, double *F, long N) {
  *U = 1.0e6;
  // Zero forces would be read as convergence: an optimiser judges a point
  // converged on force magnitude alone, so a rejected geometry with no force
  // is accepted as a minimum and the wall energy is recorded as that
  // minimum's energy. It then reaches the barrier as E_saddle - 1e6, which
  // eOn reports as a negative barrier and discards -- a real 0.062 eV copper
  // saddle was lost exactly this way. Return a force far above any
  // convergence threshold so the point can never be mistaken for a
  // stationary one, alternating its sign so the frame gains no net force.
  for (long i = 0; i < 3 * N; ++i) {
    F[i] = 0.0;
  }
  for (long i = 0; i < N; ++i) {
    F[3 * i] = (i % 2 == 0) ? 1.0 : -1.0;
  }
}

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
    // Client main arms feenableexcept unconditionally. The worker inherits
    // that mask; LAMMPS EAM (PairEAM::compute) performs IEEE divisions that
    // can raise FE_DIVBYZERO on near-coincident pairs during saddle / product
    // minimisations. Under trapping that SIGFPEs, and the continue handler
    // without MXCSR masking re-storms forever in the child (GB of identical
    // "FPE (continuing)" lines). External pot code expects soft IEEE defaults;
    // demote trapping for the whole worker process before any LAMMPS call.
    eonc::disableFPE();
    runWorkerLoop(); // never returns
  }

  // Parent: keep reqPipe write end and resPipe read end.
  //
  // Writing to a worker that has already exited raises SIGPIPE, whose default
  // action kills the client outright -- before writeExact can return the
  // error the caller is written to handle. A search that had converged on a
  // saddle at 0.062 eV died this way with signal 13 as its endpoints were
  // about to be minimised. Ignoring it turns the same condition into an
  // EPIPE return, which reaches the geometry-rejection path and respawns.
  std::signal(SIGPIPE, SIG_IGN);
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
    if (!readExact(reqFd, atomicNrs.data(),
                   sizeof(int) * static_cast<size_t>(N)) ||
        !readExact(reqFd, box, sizeof(box)) ||
        !readExact(reqFd, R.data(),
                   sizeof(double) * static_cast<size_t>(3 * N))) {
      _exit(1);
    }
    std::vector<double> mask(static_cast<size_t>(3 * N), 0.0);
    if (!readExact(reqFd, mask.data(),
                   sizeof(double) * static_cast<size_t>(3 * N))) {
      _exit(1);
    }
    setFixedMask(N, mask.data());

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
        !writeExact(resFd, F.data(),
                    sizeof(double) * static_cast<size_t>(3 * N))) {
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
    // A wedged worker never acts on the sentinel or the closed pipe, and an
    // unconditional wait then blocks the client forever at no CPU cost -- the
    // hang this teardown exists to avoid. Give the child a brief chance to
    // exit on its own, then insist.
    int st = 0;
    bool reaped = false;
    for (int i = 0; i < 100; ++i) { // up to ~1 s
      pid_t r = waitpid(workerPid, &st, WNOHANG);
      if (r == workerPid || r < 0) {
        reaped = true;
        break;
      }
      usleep(10000);
    }
    if (!reaped) {
      kill(workerPid, SIGKILL);
      waitpid(workerPid, &st, 0);
    }
    workerPid = -1;
  }
  workerSpawned = false;
}
#endif // !EONMPI && !IS_WINDOWS

void LAMMPSPot::force(long N, const double *R, const int *atomicNrs, double *F,
                      double *U, double *variance, const double *box) {
  variance = nullptr;

#ifdef EONMPI
  forceLocal(N, R, atomicNrs, F, U, box);
#elif defined(IS_WINDOWS)
  // No fork/pipe on Windows; call forceLocal directly.
  forceLocal(N, R, atomicNrs, F, U, box);
#else
  // Drive the dedicated worker process so this image's LAMMPS runs in its own
  // process (own MPI_COMM_WORLD).  Per-image NEB threads thus evaluate forces
  // as truly concurrent processes with no shared-communicator contention.
  std::lock_guard<std::mutex> workerLock(workerMutex);
  if (workerRespawnsLeft <= 0) {
    rejectGeometry(U, F, N);
    return;
  }
  ensureWorker();

  std::vector<double> mask(static_cast<size_t>(3 * N), 0.0);
  if (maskN_ == N && fixedMask_.size() == static_cast<size_t>(3 * N)) {
    mask = fixedMask_;
  }
  if (!writeExact(reqFd, &N, sizeof(N)) ||
      !writeExact(reqFd, atomicNrs, sizeof(int) * static_cast<size_t>(N)) ||
      !writeExact(reqFd, box, sizeof(double) * 9) ||
      !writeExact(reqFd, R, sizeof(double) * static_cast<size_t>(3 * N)) ||
      !writeExact(reqFd, mask.data(),
                  sizeof(double) * static_cast<size_t>(3 * N))) {
    // A worker stopped by an earlier rejected geometry leaves the request
    // pipe closed, so the first send after it fails. Respawning happens on
    // the next evaluation; reject this one rather than end the client.
    --workerRespawnsLeft;
    EONC_LOG_WARNING("[LAMMPSPot] send to worker failed; {} respawns left "
                     "(eon-7416)",
                     workerRespawnsLeft);
    stopWorker();
    rejectGeometry(U, F, N);
    return;
  }

  // eon-7416: a worker stuck on a pathological geometry (LAMMPS spinning, or
  // a NaN it never returns) would make the blocking read below hang until the
  // akmc pass times out. Bound the wait: if the worker is silent past a
  // generous per-eval deadline, kill and reap it (a stuck worker never reaches
  // EOF, so a plain waitpid would block too) and fail this evaluation so the
  // search discards the geometry; ensureWorker respawns on the next call.
  {
    struct pollfd pfd;
    pfd.fd = resFd;
    pfd.events = POLLIN;
    pfd.revents = 0;
    int pr = poll(&pfd, 1, 90000); // 90 s: orders beyond a normal force eval
    if (pr <= 0) {
      if (workerPid > 0) {
        kill(workerPid, SIGKILL);
      }
      --workerRespawnsLeft;
      EONC_LOG_WARNING(
          "[LAMMPSPot] worker force eval timed out; {} respawns left "
          "(eon-7416)",
          workerRespawnsLeft);
      stopWorker();
      rejectGeometry(U, F, N);
      return;
    }
  }
  int status = 0;
  if (!readExact(resFd, &status, sizeof(status)) ||
      !readExact(resFd, U, sizeof(double)) ||
      !readExact(resFd, F, sizeof(double) * static_cast<size_t>(3 * N))) {
    --workerRespawnsLeft;
    EONC_LOG_WARNING(
        "[LAMMPSPot] worker died during force eval; {} respawns left "
        "(eon-7416)",
        workerRespawnsLeft);
    stopWorker();
    rejectGeometry(U, F, N);
    return;
  }
  if (status != 0) {
    --workerRespawnsLeft;
    EONC_LOG_WARNING(
        "[LAMMPSPot] worker reported an evaluation error; {} respawns left "
        "(eon-7416)",
        workerRespawnsLeft);
    stopWorker();
    rejectGeometry(U, F, N);
    return;
  }
  // A saddle search that never terminates silently truncates the event
  // table: the KMC residence time is 1/sum_j k_j over the discovered
  // mechanisms, so a dropped search removes a term and biases the clock
  // (Pedersen and Jónsson, Math. Comput. Simul. 80, 1487 (2010),
  // doi:10.1016/j.matcom.2009.02.010, Fig. 1; Alexander and Schuh,
  // Modelling Simul. Mater. Sci. Eng. 24, 065014 (2016),
  // doi:10.1088/0965-0393/24/6/065014, on catalog completeness).
  // An over-aggressive saddle-search kick can drive atoms on top of
  // each other; the EAM force overflows to NaN/Inf and LAMMPS returns it
  // rather than crashing. The min-mode search then spins on non-finite
  // gradients until the akmc pass times out (0 processes). Reject the
  // evaluation so the search discards that displacement and continues.
  // Reject the geometry, not the process. Nothing between this call and
  // main catches a throw here, so the client terminates: a search that was
  // making progress is lost, and every other search sharing the pass goes
  // with it. A large finite energy with zeroed forces reads to the
  // optimiser as an impassable wall, so it backs out of the step and the
  // search abandons this configuration and carries on.
  bool nonfinite = !std::isfinite(*U);
  for (long i = 0; i < 3 * N && !nonfinite; ++i) {
    nonfinite = !std::isfinite(F[i]);
  }
  if (nonfinite) {
    EONC_LOG_WARNING("[LAMMPSPot] non-finite force or energy; rejecting "
                     "geometry (eon-7416)");
    rejectGeometry(U, F, N);
  }
#endif
}

void LAMMPSPot::forceLocal(long N, const double *R, const int *atomicNrs,
                           double *F, double *U, const double *box) {
  // Same contract as ASE / Metatomic: external pot libraries are not written
  // for FE traps. Cover in-process (EONMPI / Windows) and any path that still
  // has trapping armed when forceLocal runs. Always restore so FE traps do not
  // stay demoted for the rest of the process after the first force call.
  eonc::FPEHandler fpeh;
  fpeh.eat_fpe();
  try {
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
    applySetforce(N);
    // New instance / box change: rebuild neighbors. create_atoms sits at
    // the origin; pre no would evaluate on that neighbor list.
    if (newLammps) {
      lmp.command(LAMMPSObj, "run 1 pre yes post no");
    } else {
      lmp.command(LAMMPSObj, "run 1 pre no post no");
    }

    auto *pe =
        static_cast<double *>(lmp.extract_variable(LAMMPSObj, "pe", nullptr));
    auto *fx =
        static_cast<double *>(lmp.extract_variable(LAMMPSObj, "fx", "all"));
    auto *fy =
        static_cast<double *>(lmp.extract_variable(LAMMPSObj, "fy", "all"));
    auto *fz =
        static_cast<double *>(lmp.extract_variable(LAMMPSObj, "fz", "all"));
    if (!pe || !fx || !fy || !fz) {
      if (pe)
        free(pe);
      if (fx)
        free(fx);
      if (fy)
        free(fy);
      if (fz)
        free(fz);
      throw std::runtime_error(
          "LAMMPS: extract_variable returned null (pe/fx/fy/fz)");
    }
    *U = *pe;
    free(pe);

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
  } catch (...) {
    fpeh.restore_fpe();
    throw;
  }
  fpeh.restore_fpe();
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
  MPI_Comm_dup(mpiComm, &inst_comm); // private comm per per-image instance
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
    if (LAMMPSObj != nullptr) {
      lmp.close(LAMMPSObj);
      LAMMPSObj = nullptr;
    }
    throw std::runtime_error(
        "LAMMPS: in.lammps not found in working directory");
  }

  if (realunits) {
    lmp.command(LAMMPSObj, "units real");
  } else {
    lmp.command(LAMMPSObj, "units metal");
  }

  lmp.command(LAMMPSObj, "atom_style charge");
  lmp.command(LAMMPSObj, "atom_modify map array sort 0 0");
  lmp.command(LAMMPSObj, "neigh_modify delay 1");

  // LAMMPS restricted triclinic: (ax, by, cz, bx, cx, cy).
  // Row-major Matter cell also has ay, az, bz at box[1], box[2], box[5].
  constexpr double kTiltCut = 1.0e-8;
  if (std::abs(box[1]) > kTiltCut || std::abs(box[2]) > kTiltCut ||
      std::abs(box[5]) > kTiltCut) {
    throw std::runtime_error(
        "LAMMPS: cell is not restricted triclinic (ay/az/bz must be ~0)");
  }
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
