//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include<cassert>
#include<math.h>
#include "Morse.h"

/** @file
      @brief Morse potential for platinum
      @author Anonymous (possibly A. Pedersen or G. Henkelman ), revision: Jean Claude C. Berthet
      @date Unknown, revision: 2010, University of Iceland
      */

/** @class Morse
      @brief Morse potential for platinum.
      Default parameters are for platinum, but parameters for other atoms may be specified in setParameters()
      Morse energy function: @f$ V(r)=D_e\left[1-e^{-a(r-r_e)}\right]^2 @f$. Equation and notation from
      http://en.wikipedia.org/wiki/Morse_potential .
      */

Morse::Morse(){
      // Parameters De in eV, a in Angstroms, re in Angstroms, cutoff in Angstroms
      setParameters(0.7102, 1.6047, 2.8970, 9.5);
}

Morse::Morse(double De, double a, double re, double cutoff)
{
      setParameters(De, a, re, cutoff);
}

void Morse::cleanMemory(void){
      return;
}

void Morse::setParameters(double De, double a, double re, double cutoff)
{
      re_ = re;
      De_ = De;
      a_ = a;
      cutoff_ = cutoff;
      double f;
      morse(cutoff, energyCutoff_, f);
}

// For bindings
std::pair<AtomMatrix, double> Morse::energy_and_forces(long nAtoms, AtomMatrix positions, Matrix3d box){
      AtomMatrix forces = AtomMatrix::Constant(nAtoms, 3, 0);
      VectorXi atomicNrs = VectorXi::Constant(nAtoms, 0);
      int *atnrs = atomicNrs.data();
      double *pos = positions.data();
      double *frcs = forces.data();
      double *bx = box.data();
      double energy{0};
      Morse::force(nAtoms, pos, atnrs, frcs, &energy, bx, 1);
      return std::make_pair(forces, energy);
}

// pointer to number of atoms, pointer to array of positions	
// pointer to array of forces, pointer to internal energy
// adress to supercell size
void Morse::force(long N, const double *R, const int *atomicNrs, double *F, double *U, const double *box, int nImages=1)
{
      assert((box[0] > 0) and (box[4] > 0) and (box[8] > 0));
      double diffR=0, diffRX, diffRY, diffRZ;
      *U = 0;
      for(int i=0;i<N;i++){
            F[ 3*i ] = 0;
            F[3*i+1] = 0;
            F[3*i+2] = 0;
      };
      // Initializing en
      for(int i=0;i<N;i++){
            for(int j=i+1;j<N;j++){
                  diffRX = R[ 3*i ]-R[ 3*j ];
                  diffRY = R[3*i+1]-R[3*j+1];
                  diffRZ = R[3*i+2]-R[3*j+2];
                  diffRX = diffRX - box[0]*floor(diffRX/box[0]+0.5); // floor = largest integer value less than argument 
                  diffRY = diffRY - box[4]*floor(diffRY/box[4]+0.5);
                  diffRZ = diffRZ - box[8]*floor(diffRZ/box[8]+0.5);
                  diffR = sqrt(diffRX*diffRX+diffRY*diffRY+diffRZ*diffRZ);
                  assert(std::isnormal(diffR));
                  if(diffR < cutoff_){
                        double force, energy;
                        morse(diffR, energy, force);
                        *U+=energy;
                        F[ 3*i ]+=force*diffRX/diffR;
                        F[3*i+1]+=force*diffRY/diffR;
                        F[3*i+2]+=force*diffRZ/diffR;
                        F[ 3*j ]-=force*diffRX/diffR;
                        F[3*j+1]-=force*diffRY/diffR;
                        F[3*j+2]-=force*diffRZ/diffR;
                        *U -= energyCutoff_;
                  }
            }
      }
      return;
}

/** Calculate energy and force.
      @param[in] r Distance between two atoms
      @param[out] energy
      @param[out] force -d(energy)/dr
      @see #Morse for equation.
      */
void Morse::morse(double r, double & energy, double & force)
{
      double const d=1-exp(-a_*(r-re_));
      energy=De_*d*d-De_;
      force= 2*De_*d*(d-1)*a_;
}
