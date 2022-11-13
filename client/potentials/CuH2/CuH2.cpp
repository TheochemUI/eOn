//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include"CuH2.h"

CuH2::CuH2(Parameters *p){
}

void CuH2::initialize(void){
    return;
}

void CuH2::cleanMemory(void){
    return;
}

// pointer to number of atoms, pointer to array of positions	
// pointer to array of forces, pointer to internal energy
// address to supercell size
void CuH2::force(long N, const double *R, const int *atomicNrs, double *F, double *U, const double *box, int nImages=1){
    int natms[] {216, 2};
    int ndim {3*218};
    // The box only takes the diagonal (assumes cubic)
    double box_eam[] {box[0], box[4], box[8]};

    c_force_eam(natms, ndim, box_eam, const_cast<double*>(R), F, U);

    // for(int i=0; i<N; i++){
    //     F[ 3*i ] = fake1;
    //     F[3*i+1] = fake1;
    //     F[3*i+2] = fake1;
    // }

    // *U = fake2;
    return;
}
