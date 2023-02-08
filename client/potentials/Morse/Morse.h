//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef MORSE
#define MORSE
/** @file
      @brief Morse potential for platinum
      @author Anonymous (possibly A. Pedersen or G. Henkelman), revision: Jean Claude C. Berthet
      @date Unknown, revision: 2010, University of Iceland
      */
//#include "LJBinary.h"
#include <cmath>
#include "../../Potential.h"

class Morse : public Potential {
public:
      Morse();
      Morse(double re, double De, double a, double cutoff);
      void cleanMemory(void);// required by PotentialsInterface
      void force(long N, const double *R, const int *, double *F, double *U, const double *box);
      void initialize() {};// required by PotentialsInterface
      void setParameters(double De, double a, double re, double cutoff);
private:
      void morse(double r, double & energy, double & force);
      double re_; 
      double De_;
      double cutoff_;
      double a_;
      double energyCutoff_;
};
#endif
