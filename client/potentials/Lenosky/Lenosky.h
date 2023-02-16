//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef LENOSKY_POTENTIAL
#define LENOSKY_POTENTIAL

#include "../../Potential.h"

    /** External function implemented in Fortran; calculate interactions between atoms using Lenosky force field
    @param[in]	N           number of atoms
    @param[in]	R           array to positions of the atoms in Angstrom
    @param[out]	F           array used to return the forces between atoms, in eV/Angstrom
    @param[out]	U           pointer to energy in eV
    @param[in]  bx, by, bz  pointer to box dimensions in Angstrom
    */
extern "C" {
    void lenosky_(const long int *N, const double *R, double *F, double *U, const double* bx, const double* by, const double* bz);
}    

/** Lenosky potential */
class Lenosky : public Potential{    
public:
// Functions
	// constructor
    Lenosky(Parameters* params): Potential(params) {};
	
    // To satisfy interface
    void cleanMemory(void);
    void force(long N, const double *R, const int *atomicNrs, double *F, double *U, const double *box);
    std::pair<double, AtomMatrix> get_ef(const AtomMatrix pos,
                                         const VectorXi atmnrs,
                                         const Matrix3d box) override {
      double energy{std::numeric_limits<double>::infinity()};
      long nAtoms{pos.rows()};
      AtomMatrix forces{Eigen::MatrixXd::Zero(nAtoms, 3)};
      this->force(nAtoms, pos.data(), atmnrs.data(), forces.data(), &energy,
                  box.data());
      return std::make_pair(energy, forces);
    };
};
#endif

