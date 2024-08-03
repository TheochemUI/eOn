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

#include <cstdio>
#include <errno.h>
#include <stdlib.h>
#include <unistd.h>

#ifdef WIN32
#include <windows.h>
#define sleep(n) Sleep(1000 * n)
// #define popen _popen
#else
#include <fcntl.h>
#include <sys/wait.h>
#endif

#include "VASP.h"
namespace eonc {
bool VASP::firstRun = true;
long VASP::vaspRunCount = 0;
pid_t VASP::vaspPID = 0;

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
    dup2(outFd, 1);
    dup2(outFd, 2);
    /*
            char vaspPath[1024];
            if(strncpy(vaspPath, "vasp", 1024)) {
    //        if(strncpy("vasp", vaspPath, 1024)) {

                fprintf(stderr, "problem resolving vasp filename\n");
                exit(1);
            }
    */
    /* if (0)
            {
                if (execlp(vaspPath, "vasp", NULL) == -1)
                {
                    fprintf(stderr, "error spawning vasp: %s\n",
    strerror(errno)); exit(1);
                }
            }
            else
            {
    //            if (execlp("mpirun", "mpirun", "-n", "8", "vasp", NULL) == -1)
    */
    if (execlp("./runvasp.sh", "./runvasp.sh", NULL) == -1) {
      fprintf(stderr, "error spawning vasp: %s\n", strerror(errno));
      exit(1);
    }
    //       }
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

  // printf("vasp force call");
  // fflush(stdout);
  while (access("FU", F_OK) == -1) {
    sleep(1);
    // printf(".");
    // fflush(stdout);
    vaspRunning();
  }
  // printf("\n");
  readFU(N, efvd->F, &efvd->energy);
  remove("FU");
  vaspRunCount++;
  return;
}

void VASP::writePOSCAR(long N, const double *R, const size_t *atomicNrs,
                       const double *box) {
  // Positions are scaled
  long i = 0;
  long i_old = 0;
  FILE *POSCAR;

  POSCAR = fopen("POSCAR", "w");

  // header line (treated as a comment)
  i_old = 0;
  fprintf(POSCAR, "%d ", atomicNrs[0]);
  for (i = 0; i < N; i++) {
    if (atomicNrs[i] != atomicNrs[i_old]) {
      fprintf(POSCAR, "%d ", atomicNrs[i]);
      i_old = i;
    }
  }
  fprintf(POSCAR, ": Atomic numbers\n");

  // boundary box
  fprintf(POSCAR, "1.0\n");
  fprintf(POSCAR, " %.8f\t%.8f\t%.8f\n", box[0], box[1], box[2]);
  fprintf(POSCAR, " %.8f\t%.8f\t%.8f\n", box[3], box[4], box[5]);
  fprintf(POSCAR, " %.8f\t%.8f\t%.8f\n", box[6], box[7], box[8]);

  // the number of atoms of each different atomic type
  i_old = 0;
  for (i = 0; i < N; i++) {
    if (atomicNrs[i] != atomicNrs[i_old]) {
      fprintf(POSCAR, "%li ", i - i_old);
      i_old = i;
    }
  }
  fprintf(POSCAR, "%li\n", N - i_old);

  // coordinates for all atoms
  fprintf(POSCAR, "Cartesian\n");
  for (i = 0; i < N; i++) {
    fprintf(POSCAR, "%.19f\t%.19f\t%.19f\t T T T\n", R[i * 3 + 0], R[i * 3 + 1],
            R[i * 3 + 2]);
  }
  fclose(POSCAR);

  if (firstRun) {
    firstRun = false;
  } else {
    FILE *NEWCAR = fopen("NEWCAR", "w");
    fclose(NEWCAR);
  }

  //    FILE *NEWCAR = fopen("NEWCAR", "w");
  //    fclose(NEWCAR);

  return;
}

void VASP::readFU(long N, double *F, double *U) {
  FILE *FU;
  FU = fopen("FU", "r");

  fscanf(FU, "%lf", U);

  for (int i = 0; i < N; i++) {
    fscanf(FU, "%lf %lf %lf", &F[i * 3 + 0], &F[i * 3 + 1], &F[i * 3 + 2]);
  }
  fclose(FU);
  return;
}

} // namespace eonc
