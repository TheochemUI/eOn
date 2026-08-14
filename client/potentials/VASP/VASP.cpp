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

#include <cstddef>
#include <cstdio>
#include <cstring>
#include <errno.h>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <stdlib.h>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

#ifdef _WIN32
#include <windows.h>
#define sleep(n) Sleep(1000 * n)
// #define popen _popen
#else
#include <fcntl.h>
#include <sys/wait.h>
#include <unistd.h>
#endif

#include "eon/EonLogger.h"
#include "eon/potentials/VASP/VASP.h"

namespace {

// The driver script, and the files VASP and eOn hand back and forth. All of
// them resolve against the working directory of the calling process.
constexpr const char *kVaspScript = "runvasp.sh";
constexpr const char *kForceFile = "FU";
constexpr const char *kNewCarFile = "NEWCAR";
constexpr const char *kStopCarFile = "STOPCAR";

// The interactive handshake: eOn writes POSCAR and touches NEWCAR, VASP
// answers with FU, and eOn writes STOPCAR to end the run. Each of these left
// behind by an earlier client in the same directory breaks the next one, so
// a run clears them before starting VASP. A stale FU is read as this run's
// first forces, a stale NEWCAR feeds VASP a POSCAR eOn has not written yet,
// and a stale STOPCAR aborts VASP on its first ionic step.
constexpr const char *kHandshakeFiles[] = {kForceFile, kNewCarFile,
                                           kStopCarFile};

// Results and restart data VASP writes. VASP overwrites every one of them on
// its next run and eOn reads none of them, so removing them buys a
// calculation nothing; it discards the record of the previous one, and with
// WAVECAR, CHGCAR and TMPCAR the wavefunction and charge density a restart
// reads under ISTART and ICHARG.
constexpr const char *kStaleFiles[] = {
    "TMPCAR", "CHG",     "CHGCAR", "CONTCAR", "DOSCAR",  "EIGENVAL",
    "IBZKPT", "OSZICAR", "OUTCAR", "PCDAT",   "WAVECAR", "XDATCAR"};

// A species and how many atoms it covers.
using SpeciesRun = std::pair<int, long>;

/// \brief Split the atom list into runs of one species each.
///
/// VASP reads the species counts from one line of the POSCAR and the
/// coordinates from the block that follows, so both have to describe the same
/// ordering: every species occupies a single contiguous run. eOn does not sort
/// atoms, so an interleaved structure has to be rejected rather than written
/// out with counts that disagree with the coordinates.
std::vector<SpeciesRun> speciesRuns(long N, const int *atomicNrs) {
  std::vector<SpeciesRun> runs;
  for (long i = 0; i < N; i++) {
    if (!runs.empty() && runs.back().first == atomicNrs[i]) {
      runs.back().second++;
      continue;
    }
    for (const auto &run : runs) {
      if (run.first == atomicNrs[i]) {
        throw std::runtime_error(std::format(
            "A POSCAR needs each species in one contiguous run, but atomic "
            "number {} appears again at atom {}",
            atomicNrs[i], i));
      }
    }
    runs.emplace_back(atomicNrs[i], 1);
  }
  return runs;
}

} // namespace

bool VASP::firstRun = true;
long VASP::vaspRunCount = 0;
pid_t VASP::vaspPID = 0;

void VASP::clearHandshakeFiles() {
  for (const char *name : kHandshakeFiles) {
    std::error_code ec;
    if (std::filesystem::remove(name, ec)) {
      EONC_LOG_INFO("VASP cleared leftover {}", name);
    } else if (ec) {
      EONC_LOG_WARNING("VASP could not clear leftover {}: {}", name,
                       ec.message());
    }
  }
}

void VASP::removeStaleFiles() {
  for (const char *name : kStaleFiles) {
    std::error_code ec;
    if (std::filesystem::remove(name, ec)) {
      EONC_LOG_INFO("VASP removed leftover {}", name);
    } else if (ec) {
      EONC_LOG_WARNING("VASP could not remove leftover {}: {}", name,
                       ec.message());
    }
  }
}

void VASP::cleanMemory(void) {
  vaspRunCount--;
  if (vaspRunCount < 1) {
    // Runs from a destructor, so it reports rather than throws.
    std::ofstream stopcar(kStopCarFile, std::ios::trunc);
    if (!stopcar) {
      EONC_LOG_WARNING("Could not open {} to stop VASP", kStopCarFile);
      return;
    }
    stopcar << "LABORT = .TRUE.\n";
    stopcar.close();
    if (!stopcar) {
      EONC_LOG_WARNING("Could not write {} to stop VASP", kStopCarFile);
    }
  }
  return;
}

void VASP::spawnVASP() {
  // Runs once per client, since vaspPID stays set for the life of the
  // process. The POSCAR eOn just wrote stays, and the first call writes no
  // NEWCAR, so nothing here discards a signal this run has sent.
  clearHandshakeFiles();

  // execlp resolves a relative path against whatever directory the process
  // last changed to, which under the MPI dispatcher is not necessarily the
  // one holding the run. Resolve it in the parent, where a failure can still
  // be reported.
  std::error_code ec;
  const std::filesystem::path script =
      std::filesystem::absolute(kVaspScript, ec);
  if (ec) {
    throw std::runtime_error(
        std::format("Could not resolve {}: {}", kVaspScript, ec.message()));
  }
  if (!std::filesystem::exists(script, ec) || ec) {
    throw std::runtime_error(
        std::format("{} does not exist; VASP needs it in the working directory",
                    script.string()));
  }
  const std::string scriptPath = script.string();

  if ((vaspPID = fork()) == -1) {
    fprintf(stderr, "error forking for vasp: %s\n", strerror(errno));
    exit(1);
  }

  if (vaspPID) {
    /* We are the parent */
    setvbuf(stdout, (char *)NULL, _IONBF, 0); // non-buffered output
  } else {
    /* We are the child */
    int outFd = open("vaspout", O_CREAT | O_WRONLY | O_TRUNC, 0644);
    if (outFd == -1) {
      fprintf(stderr, "error opening vaspout: %s\n", strerror(errno));
      _exit(1);
    }
    if (dup2(outFd, 1) == -1 || dup2(outFd, 2) == -1) {
      fprintf(stderr, "error redirecting vasp output: %s\n", strerror(errno));
      _exit(1);
    }
    execl(scriptPath.c_str(), scriptPath.c_str(), (char *)NULL);
    // Only reached when the exec failed; _exit keeps the child out of the
    // parent's exit handlers.
    fprintf(stderr, "error spawning vasp: %s\n", strerror(errno));
    _exit(1);
  }
}

bool VASP::vaspRunning() {
  pid_t pid;
  int status;

  if (vaspPID == 0) {
    return false;
  }

  pid = waitpid(vaspPID, &status, WNOHANG);

  if (pid) {
    fprintf(stderr, "vasp died unexpectedly!\n");
    exit(1);
  }

  return true;
}

void VASP::force(long N, const double *R, const int *atomicNrs, double *F,
                 double *U, double *variance, const double *box) {
  variance = nullptr;
  writePOSCAR(N, R, atomicNrs, box);

  if (!vaspRunning()) {
    spawnVASP();
  }

  // access() succeeds the moment VASP creates FU, not when it finishes
  // writing it, so this can observe a partial file. readFU treats a short
  // read as an error rather than absorbing it, which turns the race into a
  // failure instead of wrong forces. Closing it properly needs a completion
  // signal from runvasp.sh, such as writing to a temporary name and renaming
  // it once the write completes.
  while (access(kForceFile, F_OK) == -1) {
    sleep(1);
    vaspRunning();
  }
  readFU(N, F, U);

  std::error_code ec;
  if (!std::filesystem::remove(kForceFile, ec) || ec) {
    // Leaving it behind means the next call reads this call's numbers.
    throw std::runtime_error(std::format("Could not remove {}: {}", kForceFile,
                                         ec ? ec.message() : "already gone"));
  }
  vaspRunCount++;
  return;
}

void VASP::writePOSCAR(long N, const double *R, const int *atomicNrs,
                       const double *box) {
  // Positions are scaled
  const std::vector<SpeciesRun> runs = speciesRuns(N, atomicNrs);

  std::ofstream poscar("POSCAR", std::ios::trunc);
  if (!poscar) {
    throw std::runtime_error("Could not open POSCAR for writing");
  }

  // header line (treated as a comment)
  for (const auto &run : runs) {
    poscar << std::format("{} ", run.first);
  }
  poscar << ": Atomic numbers\n";

  // boundary box
  poscar << "1.0\n";
  for (int i = 0; i < 3; i++) {
    poscar << std::format(" {:.8f}\t{:.8f}\t{:.8f}\n", box[i * 3 + 0],
                          box[i * 3 + 1], box[i * 3 + 2]);
  }

  // the number of atoms of each different atomic type
  for (std::size_t i = 0; i < runs.size(); i++) {
    poscar << std::format("{}{}", runs[i].second,
                          i + 1 == runs.size() ? "\n" : " ");
  }

  // coordinates for all atoms
  poscar << "Cartesian\n";
  for (long i = 0; i < N; i++) {
    poscar << std::format("{:.19f}\t{:.19f}\t{:.19f}\t T T T\n", R[i * 3 + 0],
                          R[i * 3 + 1], R[i * 3 + 2]);
  }

  poscar.close();
  if (!poscar) {
    throw std::runtime_error("Could not write the structure to POSCAR");
  }

  if (firstRun) {
    firstRun = false;
  } else {
    // An empty NEWCAR tells the running VASP that a new POSCAR is ready.
    std::ofstream newcar(kNewCarFile, std::ios::trunc);
    if (!newcar) {
      throw std::runtime_error(
          std::format("Could not open {} to signal VASP", kNewCarFile));
    }
    newcar.close();
    if (!newcar) {
      throw std::runtime_error(
          std::format("Could not write {} to signal VASP", kNewCarFile));
    }
  }

  return;
}

void VASP::readFU(long N, double *F, double *U) {
  std::ifstream fu(kForceFile);
  if (!fu) {
    throw std::runtime_error(
        std::format("Could not open {}; VASP left no result", kForceFile));
  }

  if (!(fu >> *U)) {
    throw std::runtime_error(
        std::format("Could not read the energy from {}", kForceFile));
  }

  for (long i = 0; i < N; i++) {
    if (!(fu >> F[i * 3 + 0] >> F[i * 3 + 1] >> F[i * 3 + 2])) {
      throw std::runtime_error(std::format(
          "{} holds forces for {} atoms, expected {}", kForceFile, i, N));
    }
  }
  return;
}
