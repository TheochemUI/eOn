//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include "GPAW.h"
#include <mpi.h>
#include <stdio.h>
#include <unistd.h>

int GPAW::instances = 0;
bool GPAW::firstRun = true;
MPI::Intercomm GPAW::gpawComm = MPI::COMM_NULL;

GPAW::GPAW(void)
{
    if (firstRun) {
        MPI::Init();
        firstRun = false;
    }

    if (instances == 0) {
        //XXX: hard coded number of processes
        int numProcs = 2;
        gpawComm = MPI::COMM_WORLD.Spawn("gpaw_sp.py",
                       MPI::ARGV_NULL, numProcs,
                       MPI::INFO_NULL, 0,
                       MPI_ERRCODES_IGNORE);
    }
    instances++;
    return;
}

void GPAW::cleanMemory(void)
{
    instances--;
    if (instances == 0) {
        int die=1;
        gpawComm.Send(&die, 1, MPI::INT, 0, 0);
        gpawComm.Disconnect();
    }
    return;
}

GPAW::~GPAW()
{
	cleanMemory();
}

void GPAW::force(long N, const double *R, const int *atomicNrs, double *F, 
                 double *U, const double *box)
{
    int die=0;
    gpawComm.Send(&die, 1, MPI::INT, 0, 0);
    gpawComm.Send(&N, 1, MPI::LONG, 0, 0);
    gpawComm.Send(atomicNrs, N, MPI::INT, 0, 0);
    gpawComm.Send(R, 3*N, MPI::DOUBLE, 0, 0);
    gpawComm.Send(box, 9, MPI::DOUBLE, 0, 0);
    int pbc=1;
    gpawComm.Send(&pbc, 1, MPI::INT, 0, 0);
    MPI::Status status;
    int failed;
    gpawComm.Recv(&failed, 1, MPI::INT, 0, 0, status);
    if (failed == 1) {
        throw 123;
    }
    gpawComm.Recv(U, 1, MPI::DOUBLE, 0, 0, status);
    gpawComm.Recv(F, 3*N, MPI::DOUBLE, 0, 0, status);
}
