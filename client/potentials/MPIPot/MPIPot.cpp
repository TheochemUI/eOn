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

#include "MPIPot.h"
#include <mpi.h>
#include <stdio.h>
#include <time.h>
#ifndef _WIN32
#include <unistd.h>
#endif

MPIPot::MPIPot(const Parameters &p)
    : Potential(p) {
  potentialRank = p.potential_options.MPIPotentialRank;
  poll_period = p.potential_options.MPIPollPeriod;
  return;
}

void MPIPot::cleanMemory(void) { return; }

MPIPot::~MPIPot() { cleanMemory(); }

void MPIPot::force(long N, const double *R, const int *atomicNrs, double *F,
                   double *U, double *variance, const double *box) {
  variance = nullptr;
  // Send data to potential
  int pbc = 1;
  int failed;
  char cwd[1024];
  long icwd[1024];
  getcwd(cwd, 1024);
  for (int i = 0; i < 1024; i++) {
    icwd[i] = static_cast<long>(cwd[i]);
  }
  int intn = static_cast<int>(N);
  MPI_Send(&intn, 1, MPI_INT, potentialRank, 0, MPI_COMM_WORLD);
  MPI_Send(atomicNrs, N, MPI_INT, potentialRank, 0, MPI_COMM_WORLD);
  MPI_Send(R, 3 * N, MPI_DOUBLE, potentialRank, 0, MPI_COMM_WORLD);
  MPI_Send(box, 9, MPI_DOUBLE, potentialRank, 0, MPI_COMM_WORLD);
  MPI_Send(&pbc, 1, MPI_INT, potentialRank, 0, MPI_COMM_WORLD);
  MPI_Send(&icwd[0], 1024, MPI_INT, potentialRank, 0, MPI_COMM_WORLD);

  if (poll_period > 0.0) {
    int eon_flag = 0;
    MPI_Iprobe(potentialRank, 0, MPI_COMM_WORLD, &eon_flag, MPI_STATUS_IGNORE);
    while (!eon_flag) {
      usleep(static_cast<useconds_t>(poll_period / 1000000.0));
      MPI_Iprobe(potentialRank, 0, MPI_COMM_WORLD, &eon_flag, MPI_STATUS_IGNORE);
    }
  }

  // Recv data from potential
  MPI_Recv(&failed, 1, MPI_INT, potentialRank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  if (failed == 1) {
    throw 100;
  }

  MPI_Recv(U, 1, MPI_DOUBLE, potentialRank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  MPI_Recv(F, 3 * N, MPI_DOUBLE, potentialRank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  // printf("energy: %12.4e\n", *U);
  // printf("forces:\n");
  // for (int i=0;i<N;i++) printf("%12.4e %12.4e %12.4e\n", F[3*i], F[3*i+1],
  // F[3*i+2]);
}
