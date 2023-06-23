//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

// serves as an interface between emt potentials provided by CamposASE and
// dynamics provided by EON

#ifndef EFFECTIVE_MEDIUM_THEORY
#define EFFECTIVE_MEDIUM_THEORY

#include "Asap/Atoms.h"
#include "Asap/EMT.h"
#include "Asap/EMTDefaultParameterProvider.h"
#include "Asap/EMTRasmussenParameterProvider.h"
#include "Asap/SuperCell.h"
#include "Asap/Vec.h"

#include "../../Parameters.h"
#include "../../Potential.h"

/** EMT potential. Inspect the EMT_parms.h to see what the EMT potential is
 * hardcoded to describe.*/
class EffectiveMediumTheory : public Potential {

private:
  //	Variables
  long numberOfAtoms;
  bool periodicity[3];
  Atoms *AtomsObj;
  EMTDefaultParameterProvider *EMTParameterObj;
  EMT *EMTObj;
  SuperCell *SuperCellObj;

public:
  // Functions
  // constructor and destructor
  EffectiveMediumTheory(std::shared_ptr<Parameters> p) : Potential(p) {
    // dummy variables
    AtomsObj = 0;
    EMTObj = 0;
    SuperCellObj = 0;
    EMTParameterObj = 0;
    numberOfAtoms = 0;

    // should have periodic boundary conditions in all directions
    periodicity[0] = true;
    periodicity[1] = true;
    periodicity[2] = true;
  };
  ~EffectiveMediumTheory(void){};
  void cleanMemory(void);

  // To satify interface
  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box);
};

#endif
