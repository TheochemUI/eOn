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

#include "VASP.h"

#include <cerrno>
#include <cstring>
#include <fcntl.h>
#include <iostream>
#include <stdexcept>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

namespace eonc {

bool VASP::firstRun = true;
long VASP::vaspRunCount = 0;
pid_t VASP::vaspPID = 0;

void runVaspProcess() {
  int outFd = open("vaspout", O_CREAT | O_WRONLY | O_TRUNC, 0644);
  if (outFd == -1) {
    throw std::runtime_error(std::string("error opening vaspout: ") +
                             std::strerror(errno));
  }
  dup2(outFd, STDOUT_FILENO);
  dup2(outFd, STDERR_FILENO);
  close(outFd);

  if (execlp("./runvasp.sh", "./runvasp.sh", nullptr) == -1) {
    throw std::runtime_error(std::string("error spawning vasp: ") +
                             std::strerror(errno));
  }
}

void VASP::cleanMemory(void) {
  vaspRunCount--;
  if (vaspRunCount < 1) {
    FILE *stopcar = fopen("STOPCAR", "w");
    fprintf(stopcar, "LABORT = .TRUE.\n");
    fclose(stopcar);
  }
  return;
}

void VASP::spawnVASP() {
  vaspProcess = std::async(std::launch::async, []() {
    try {
      runVaspProcess();
    } catch (const std::exception &e) {
      std::cerr << e.what() << std::endl;
      std::exit(EXIT_FAILURE);
    }
  });
}

bool VASP::vaspRunning() {
  if (vaspPID == 0) {
    return false;
  }
  int status;
  pid_t pid = waitpid(vaspPID, &status, WNOHANG);

  if (pid == -1) {
    std::cerr << "Error in waitpid\n";
    std::exit(EXIT_FAILURE);
  } else if (pid > 0) {
    std::cerr << "vasp died unexpectedly!\n";
    std::exit(EXIT_FAILURE);
  }
  return true;
}

void VASP::forceImpl(const ForceInput &fip, ForceOut *efvd) {
#ifdef EON_CHECKS
  eonc::pot::checkParams(fip);
  eonc::pot::zeroForceOut(fip.nAtoms, efvd);
#endif
  const long int N = fip.nAtoms;
  writePOSCAR(N, fip.pos, fip.atmnrs, fip.box);

  if (!vaspRunning()) {
    spawnVASP();
  }

  // Wait for the "FU" file to be created
  while (!std::filesystem::exists("FU")) {
    std::this_thread::sleep_for(std::chrono::milliseconds(1));
    vaspRunning();
  }

  readFU(N, efvd->F, &efvd->energy);
  std::filesystem::remove("FU");
  vaspRunCount++;
}

void VASP::writePOSCAR(long N, const double *R, const size_t *atomicNrs,
                       const double *box) {
  // Positions are scaled
  long i = 0;
  long i_old = 0;

  std::ofstream poscar("POSCAR");

  if (!poscar.is_open()) {
    throw std::runtime_error("Could not open POSCAR file for writing.");
  }

  // header line (treated as a comment)
  std::string header = fmt::format("{} ", atomicNrs[0]);
  i_old = 0;
  for (i = 0; i < N; i++) {
    if (atomicNrs[i] != atomicNrs[i_old]) {
      header += fmt::format("{} ", atomicNrs[i]);
      i_old = i;
    }
  }
  header += ": Atomic numbers\n";
  poscar << header;

  // boundary box
  poscar << "1.0\n";
  for (i = 0; i < 9; i += 3) {
    poscar << fmt::format(" {:.8f}\t{:.8f}\t{:.8f}\n", box[i], box[i + 1],
                          box[i + 2]);
  }

  // the number of atoms of each different atomic type
  i_old = 0;
  for (i = 0; i < N; i++) {
    if (atomicNrs[i] != atomicNrs[i_old]) {
      poscar << fmt::format("{} ", i - i_old);
      i_old = i;
    }
  }
  poscar << fmt::format("{}\n", N - i_old);

  // coordinates for all atoms
  poscar << "Cartesian\n";
  for (i = 0; i < N; i++) {
    poscar << fmt::format("{:.19f}\t{:.19f}\t{:.19f}\t T T T\n", R[i * 3 + 0],
                          R[i * 3 + 1], R[i * 3 + 2]);
  }
  poscar.close();

  if (firstRun) {
    firstRun = false;
  } else {
    std::ofstream newcar("NEWCAR");
    if (!newcar.is_open()) {
      throw std::runtime_error("Could not open NEWCAR file for writing.");
    }
    newcar.close();
  }
}

void VASP::readFU(long N, double *F, double *U) {
  std::ifstream fu("FU");

  if (!fu.is_open()) {
    throw std::runtime_error("Could not open FU file for reading.");
  }

  fu >> *U;

  for (long i = 0; i < N; ++i) {
    fu >> F[i * 3 + 0] >> F[i * 3 + 1] >> F[i * 3 + 2];
  }

  if (!fu) {
    throw std::runtime_error("Error reading FU file.");
  }

  fu.close();
}

} // namespace eonc
