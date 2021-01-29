//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include"GPRPotential.h"
#include"../../../gpr_dimer/structures/Structures.h"
#include"../../../gpr_dimer/gpr/auxiliary/AdditionalFunctionality.h"

GPRPotential::GPRPotential(Parameters *p){
    gpr_model = nullptr;
}

void GPRPotential::registerGPRObject(GaussianProcessRegression *_gpr_model){
    gpr_model = _gpr_model;
}

void GPRPotential::initialize(void){
}

void GPRPotential::cleanMemory(void){
}

// pointer to number of atoms, pointer to array of positions	
// pointer to array of forces, pointer to internal energy
// adress to supercell size
void GPRPotential::force(long N, const double *R, const int *atomicNrs, double *F, double *U, const double *box){
    Observation observation;

    // Copy R points. Note, R should correspond to the moving atoms only.
    observation.R.resize(1, N * 3);
    for(int i=0; i<N; i++){
        observation.R.set(i, {R[ 3*i ], R[3*i+1], R[3*i+2]});
    }

    // Note, the following functions should be called before calling for gpr_model->calculatePotential()
    // gpr_model->decomposeCovarianceMatrix(R, ind) - takes covariance matrix and vector of repetitive indices
    // gpr_model->calculateMeanPrediction() - takes a vector of combined energy and force
    // gpr_model->calculatePosteriorMeanPrediction() - no arguments
    gpr_model->calculatePotential(observation);

    for(int i=0; i<N; i++){
        F[ 3*i ] = observation.G[ 3*i ];
        F[3*i+1] = observation.G[3*i+1];
        F[3*i+2] = observation.G[3*i+2];
    }
    
    *U = observation.E;
}
