// serves as an interface between LAMMPS potentials maintained by SANDIA

#pragma once

#include "../../Parameters.h"
#include "../../Potential.h"

class LAMMPSPot : public Potential {

public:
  LAMMPSPot(std::shared_ptr<Parameters> p);
  ~LAMMPSPot(void);
  void initialize() {};
  void cleanMemory(void);
  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box);

private:
  long numberOfAtoms;
  double oldBox[9];
  void *LAMMPSObj;
  void makeNewLAMMPS(long N, const double *R, const int *atomicNrs,
                     const double *box);
  bool realunits;
};
