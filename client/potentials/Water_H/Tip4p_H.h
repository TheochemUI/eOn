//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef TIP4_H_POTENTIAL
#define TIP4_H_POTENTIAL

#include "../../Potential.h"
#include "../Water/Water.hpp"

/** External function implemented in Fortran. Calculate interactions between
atoms using the potH-H2O force field.
@param[in]	R           Array of length of the three vectors connecting the
H atom to the individual atoms in the water molecule in Angstrom. R(0)=H-O,
R(1)=H-H1 and R(2)=H-H2
@param[out]	U           Pointer to energy in kJ/mol.
@param[out]	F           Array of the norms for the three forces in
kJ/mol/Angstrom. F(0)=H-O, F(1)=H-H1 and F(2)=H-H2. Forces are . Direction of
the vectors for the distances and forces are equvalent
*/
extern "C" {
void poth2oh_(double *R, double *U, double *F);
void setup_(void);
}

/** H-Water potential.*/
class Tip4p_H : public Potential {

public:
  // Functions
  // constructor
  Tip4p_H(std::shared_ptr<Parameters> params) : Potential(params) {
    setup_();
    tip4p_pot = std::make_shared<Tip4p>(params);
  };

  // To satisfy interface
  void initialize(void);
  void cleanMemory(void);
  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box);

private:
  std::shared_ptr<Tip4p> tip4p_pot;
};
#endif
