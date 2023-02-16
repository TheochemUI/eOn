//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef NEWPOT_INTERFACE
#define NEWPOT_INTERFACE

#include "../../Potential.h"

/** Template to use if user want to provide potential. */
class NewPot : public Potential {

private:
  //	Variables
  double fake1;
  double fake2;

public:
  // Functions
  // constructor and destructor
  NewPot(Parameters *p): Potential(p), fake1{0}, fake2{0} {};

  // To satisfy interface
  void cleanMemory(void);

  std::pair<double, AtomMatrix> get_ef(const AtomMatrix pos,
                                       const VectorXi atmnrs,
                                       const Matrix3d m_box) override;
};
#endif
