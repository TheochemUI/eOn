//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include "MPIPot.h"
#include <mpi.h>
#include <stdio.h>
#include <time.h>
#include <unistd.h>

MPIPot::MPIPot(std::shared_ptr<Parameters> p) : Potential(p) {
  potentialRank = p->MPIPotentialRank;
  poll_period = p->MPIPollPeriod;
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
    icwd[i] = (long)cwd[i];
  }
  int intn = (int)N;
  MPI::COMM_WORLD.Send(&intn, 1, MPI::INT, potentialRank, 0);
  MPI::COMM_WORLD.Send(atomicNrs, N, MPI::INT, potentialRank, 0);
  MPI::COMM_WORLD.Send(R, 3 * N, MPI::DOUBLE, potentialRank, 0);
  MPI::COMM_WORLD.Send(box, 9, MPI::DOUBLE, potentialRank, 0);
  MPI::COMM_WORLD.Send(&pbc, 1, MPI::INT, potentialRank, 0);
  MPI::COMM_WORLD.Send(&icwd[0], 1024, MPI::INT, potentialRank, 0);

  if (poll_period > 0.0) {
    while (MPI::COMM_WORLD.Iprobe(potentialRank, 0) == false) {
      usleep((useconds_t)(poll_period / 1000000.0));
    }
  }

  // Recv data from potential
  MPI::COMM_WORLD.Recv(&failed, 1, MPI::INT, potentialRank, 0);
  if (failed == 1) {
    throw 100;
  }

  MPI::COMM_WORLD.Recv(U, 1, MPI::DOUBLE, potentialRank, 0);
  MPI::COMM_WORLD.Recv(F, 3 * N, MPI::DOUBLE, potentialRank, 0);
  // printf("energy: %12.4e\n", *U);
  // printf("forces:\n");
  // for (int i=0;i<N;i++) printf("%12.4e %12.4e %12.4e\n", F[3*i], F[3*i+1],
  // F[3*i+2]);
}
