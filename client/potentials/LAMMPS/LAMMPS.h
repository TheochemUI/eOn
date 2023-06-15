// serves as an interface between LAMMPS potentials maintained by SANDIA

#ifndef LAMMPS
#define LAMMPS

#include "../../Parameters.h"
#include "../../Potential.h"

class lammps : public Potential {

public:
  lammps(std::shared_ptr<Parameters> p);
  ~lammps(void);
  void initialize(){};
  void cleanMemory(void);
  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, const double *box);

private:
  long numberOfAtoms;
  double oldBox[9];
  void *LAMMPSObj;
  void makeNewLAMMPS(long N, const double *R, const int *atomicNrs,
                     const double *box);
  bool realunits;
};
#endif
