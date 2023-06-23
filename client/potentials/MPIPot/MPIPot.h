//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef __MPI_POTENTIAL__
#define __MPI_POTENTIAL__

#include "../../Parameters.h"
#include "../../Potential.h"

class MPIPot : public Potential {

public:
  MPIPot(std::shared_ptr<Parameters> p);
  ~MPIPot();
  void initialize(){};
  void cleanMemory(void);
  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box);

private:
  int potentialRank;
  double poll_period;
};

#endif
