//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef LENNARD_JONES
#define LENNARD_JONES

#include <iostream>
#include <math.h>
// #include "../../system_unit.h" // unit converters
#include "../../Potential.h"

/** Lennard Jones potential.*/
class LJ : public Potential {
private:
  double u0;
  double cuttOffR;
  double psi;
  double cuttOffU;

public:
  LJ(std::shared_ptr<Parameters> params)
      : LJ(PotType::LJ, params, 1.0, 15.0, 1.0) {}

  LJ(PotType ptype, std::shared_ptr<Parameters> params, double u0,
     double cuttOffR, double psi)
      : Potential(ptype, params),
        u0{u0},
        cuttOffR{cuttOffR},
        psi{psi} {}

  ~LJ() = default;

  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box) override;

  void setParameters(double r0Recieved, double u0Recieved, double psiRecieved);
};
#endif
