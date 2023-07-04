//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef VASP_POTENTIAL
#define VASP_POTENTIAL

#include "../../Potential.h"

class VASP : public Potential {

public:
  VASP(shared_ptr<Parameters> p) : Potential(p) {
    vaspRunCount++;
    // deleting leftovers from previous run
    system("rm -f TMPCAR");
    system("rm -f CHG");
    system("rm -f CHGCAR");
    system("rm -f CONTCAR");
    system("rm -f DOSCAR");
    system("rm -f EIGENVAL");
    system("rm -f IBZKPT");
    system("rm -f NEWCAR");
    system("rm -f FU");
    system("rm -f OSZICAR");
    system("rm -f OUTCAR");
    system("rm -f PCDAT");
    system("rm -f POSCAR");
    system("rm -f TMPCAR");
    system("rm -f WAVECAR");
    system("rm -f XDATCAR");
  }
  ~VASP() { cleanMemory(); }
  void initialize(){};
  void cleanMemory(void);
  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box);

private:
  void writePOSCAR(long N, const double *R, const int *atomicNrs,
                   const double *box);
  void readFU(long N, double *F, double *U);
  void spawnVASP();
  bool vaspRunning();
  static bool firstRun;
  static long vaspRunCount;
  static pid_t vaspPID;
};

#endif
