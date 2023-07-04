//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef AMS_IO_POT
#define AMS_IO_POT

#include "../../Matter.h"
#include "../../Potential.h"

class AMS_IO : public Potential {

public:
  AMS_IO(std::shared_ptr<Parameters> p);
  ~AMS_IO();
  void initialize(){};
  void cleanMemory(void);
  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box);

private:
  void passToSystem(long N, const double *R, const int *atomicNrs,
                    const double *box);
  void recieveFromSystem(long N, double *F, double *U);
  const char *engine;
  const char *model;
  const char *forcefield;
  const char *xc;
};

#endif
