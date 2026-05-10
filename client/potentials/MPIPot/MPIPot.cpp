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

#ifndef _WIN32
#include <unistd.h>
#endif

// MPI C++ bindings (MPI::COMM_WORLD, MPI::INT, MPI::DOUBLE) were
// deprecated in MPI-2.2 and removed in MPI-3.0 (~2012); modern
// conda-forge MPICH 4.x and OpenMPI 5.x ship without mpicxx.h. The
// translations below use the C bindings, which every MPI
// implementation guarantees:
//   MPI::COMM_WORLD.Send(...)   ->  MPI_Send(..., MPI_COMM_WORLD)
//   MPI::COMM_WORLD.Recv(...)   ->  MPI_Recv(..., MPI_COMM_WORLD,
//                                            MPI_STATUS_IGNORE)
//   MPI::COMM_WORLD.Iprobe(...) ->  MPI_Iprobe(..., MPI_COMM_WORLD,
//                                              &flag, MPI_STATUS_IGNORE)
//   MPI::INT / MPI::DOUBLE      ->  MPI_INT / MPI_DOUBLE
//
// This is a pure compile fix; it doesn't change wire semantics.
// Runtime-loadable libmpi via an MpiLoader (the FlexiBLAS-of-MPI
// pattern) is a separate follow-up; see the ABI note in
// MPIPot's docs/source/user_guide/mpi_potential.md.

MPIPot::MPIPot(const Parameters &p)
    : Potential(p) {
  potentialRank = p.potential_options.MPIPotentialRank;
  poll_period = p.potential_options.MPIPollPeriod;
}

void MPIPot::cleanMemory(void) {}

MPIPot::~MPIPot() { cleanMemory(); }

void MPIPot::force(long N, const double *R, const int *atomicNrs, double *F,
                   double *U, double *variance, const double *box) {
  variance = nullptr;
  // Send data to potential
  int pbc = 1;
  int failed = 0;
  char cwd[1024];
  // Wire format pre-dates this commit: 1024 MPI_INTs holding char
  // values (cwd[i] cast through int). The original C++-bindings code
  // declared a long[1024] but tagged the message MPI::INT, which on
  // little-endian sent the low 4 bytes of each long. Switching the
  // local buffer to int[] keeps the wire format identical on every
  // endianness (any deployed MPI server still parses it correctly).
  int icwd[1024];
  if (getcwd(cwd, sizeof(cwd)) == nullptr) {
    cwd[0] = '\0';
  }
  for (int i = 0; i < 1024; i++) {
    icwd[i] = static_cast<int>(cwd[i]);
  }
  int intn = static_cast<int>(N);
  MPI_Send(&intn, 1, MPI_INT, potentialRank, 0, MPI_COMM_WORLD);
  MPI_Send(const_cast<int *>(atomicNrs), static_cast<int>(N), MPI_INT,
           potentialRank, 0, MPI_COMM_WORLD);
  MPI_Send(const_cast<double *>(R), static_cast<int>(3 * N), MPI_DOUBLE,
           potentialRank, 0, MPI_COMM_WORLD);
  MPI_Send(const_cast<double *>(box), 9, MPI_DOUBLE, potentialRank, 0,
           MPI_COMM_WORLD);
  MPI_Send(&pbc, 1, MPI_INT, potentialRank, 0, MPI_COMM_WORLD);
  MPI_Send(icwd, 1024, MPI_INT, potentialRank, 0, MPI_COMM_WORLD);

  if (poll_period > 0.0) {
    int flag = 0;
    do {
      MPI_Iprobe(potentialRank, 0, MPI_COMM_WORLD, &flag, MPI_STATUS_IGNORE);
      if (!flag) {
        usleep(static_cast<useconds_t>(poll_period / 1000000.0));
      }
    } while (!flag);
  }

  // Recv data from potential
  MPI_Recv(&failed, 1, MPI_INT, potentialRank, 0, MPI_COMM_WORLD,
           MPI_STATUS_IGNORE);
  if (failed == 1) {
    throw 100;
  }

  MPI_Recv(U, 1, MPI_DOUBLE, potentialRank, 0, MPI_COMM_WORLD,
           MPI_STATUS_IGNORE);
  MPI_Recv(F, static_cast<int>(3 * N), MPI_DOUBLE, potentialRank, 0,
           MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}
