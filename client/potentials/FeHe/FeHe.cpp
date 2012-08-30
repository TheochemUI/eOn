//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//
//-----------------------------------------------------------------------------------
#include "FeHe.h"
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

FeHe::FeHe(void)
{
    return;
}

void FeHe::initialize(void)
{
//    potinit_();
    return;
}

FeHe::~FeHe()
{
}

// pointer to number of atoms, pointer to array of positions
// pointer to array of forces, pointer to internal energy
// address to supercell size
void FeHe::force(long N, const double *R, const int *atomicNrs, double *F, double *U, const double *box)
{
    assert((int)N > 1);

    long i;
    double *RX = new double[N];
    double *RY = new double[N];
    double *RZ = new double[N];
    double *FX = new double[N];
    double *FY = new double[N];
    double *FZ = new double[N];
    long *ISPEC = new long[N];

    for(i=0; i<N; i++){
        RX[i] = R[i*N];
        RY[i] = R[i*N+1];
        RZ[i] = R[i*N+2];
        assert(atomicNrs[i]==26 || atomicNrs[i]==2);
        if(atomicNrs[i] == 26){
            ISPEC[i] = 0;
        }else{
            ISPEC[i] = 1;
        }
    }

    feforce_(&N, RX, RY, RZ, ISPEC, FX, FY, FZ, U, &box[0], &box[4], &box[8]);

    for(i=0; i<N; i++){
        F[i*N]   = FX[i];
        F[i*N+1] = FY[i];
        F[i*N+2] = FZ[i];
    }

    delete [] RX;
    delete [] RY;
    delete [] RZ;
    delete [] FX;
    delete [] FY;
    delete [] FZ;
    delete [] ISPEC;
    return;
}
