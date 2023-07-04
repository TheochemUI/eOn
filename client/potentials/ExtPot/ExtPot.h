//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef EXT_POT
#define EXT_POT

#include "../../Potential.h"

class ExtPot : public Potential {

public:
  ExtPot(std::shared_ptr<Parameters> p)
      : Potential(p), eon_extpot_path{p->extPotPath.c_str()} {};
  ~ExtPot();
  void cleanMemory(void);
  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box) override;

private:
  void passToSystem(long N, const double *R, const int *atomicNrs,
                    const double *box);
  void recieveFromSystem(long N, double *F, double *U);
  const char *eon_extpot_path;
};

#endif
