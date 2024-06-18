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
#ifndef ALUMINUM_POTENTIAL
#define ALUMINUM_POTENTIAL

#include "../../Potential.h"

/** External function implemented in Fortran
@param[in]	N           number of atoms
@param[in]	R           array to positions of the atoms in Angstrom
@param[out]	F           array used to return the forces between atoms, in
eV/Angstrom
@param[out]	U           pointer to energy in eV
@param[in]  bx, by, bz  pointer to box dimensions in Angstrom
*/

extern "C" {
void force_(const long int *N, const double *R, double *F, double *U,
            const double *bx, const double *by, const double *bz);
}
extern "C" {
void potinit_();
}

/** Aluminum potential.*/
class Aluminum : public Potential {
public:
  Aluminum(std::shared_ptr<Parameters> params)
      : Potential(PotType::EAM_AL, params) {
    potinit_();
  };
  ~Aluminum(void){};
  // To satisfy interface
  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box) override;
};
#endif
