//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef IMD_POTENTIAL
#define IMD_POTENTIAL

#include "../../Potential.h"

class IMD : public Potential {

public:
  IMD(std::shared_ptr<Parameters> params) : Potential(params){};
  ~IMD();
  void cleanMemory(void);
  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box) override;

private:
  void writeConfIMD(long N, const double *R, const int *atomicNrs,
                    const double *box);
  void readForceIMD(long N, double *F, double *U);
  //        void spawnVASP();
  //        bool vaspRunning();
  //        static bool firstRun;
  //        static long vaspRunCount;
  //        static pid_t vaspPID;
};

#endif
