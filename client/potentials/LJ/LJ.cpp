//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include"LJ.h"

LJ::LJ(){
    this->setParameters(1.0, 15.0, 1.0);
}

LJ::LJ(double u0Recieved, double cuttOffRRecieved, double psiRecieved){
    this->setParameters(u0Recieved, cuttOffRRecieved, psiRecieved);
    return;
}

void LJ::cleanMemory(void){
    return;
}

// General Functions
void LJ::setParameters(double u0Recieved, double cuttOffRRecieved, double psiRecieved){
    u0 = u0Recieved;
    psi = psiRecieved;
    
    cuttOffR = cuttOffRRecieved;
    cuttOffU = 4*u0*(pow(psi/cuttOffR,12)-pow(psi/cuttOffR,6));
    return;
}

// pointer to number of atoms, pointer to array of positions	
// pointer to array of forces, pointer to internal energy
// adress to supercell size
void LJ::force(long N, const double *R, const int *atomicNrs, double *F, double *U, const double *box){
    double diffR=0, diffRX, diffRY, diffRZ, dU, a, b;
    *U = 0;    
    for(int i=0;i<N;i++){
        F[ 3*i ] = 0;
        F[3*i+1] = 0;
        F[3*i+2] = 0;
    }
    // Initializing end
    
    for(int i=0; i<N-1; i++){
        for(int j=i+1; j<N ;j++){
            diffRX = R[ 3*i ]-R[ 3*j ];
            diffRY = R[3*i+1]-R[3*j+1];
            diffRZ = R[3*i+2]-R[3*j+2];

            diffRX = diffRX-box[0]*floor(diffRX/box[0]+0.5); // floor = largest integer value less than argument 
            diffRY = diffRY-box[4]*floor(diffRY/box[4]+0.5);
            diffRZ = diffRZ-box[8]*floor(diffRZ/box[8]+0.5);
            
            diffR = sqrt(diffRX*diffRX+diffRY*diffRY+diffRZ*diffRZ);
                
            if(diffR<cuttOffR){                
                // 4u0((psi/r0)^12-(psi/r0)^6)
                a = pow(psi/diffR,6);
                b = 4*u0*a;
                
                *U = *U+b*(a-1)-cuttOffU;

                dU=-6*b/diffR*(2*a-1);
                // F is the negative derivative
                F[ 3*i ]=F[ 3*i ] - dU*diffRX/diffR;
                F[3*i+1]=F[3*i+1] - dU*diffRY/diffR;
                F[3*i+2]=F[3*i+2] - dU*diffRZ/diffR;
                
                F[ 3*j ]=F[ 3*j ] + dU*diffRX/diffR;
                F[3*j+1]=F[3*j+1] + dU*diffRY/diffR;
                F[3*j+2]=F[3*j+2] + dU*diffRZ/diffR;
            }
        }
    }
    return;
}

LJ::~LJ()
{
    cleanMemory();
}
