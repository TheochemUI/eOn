#ifndef MORSE
#define MORSE
/** @file
      @brief Morse potential for platinum
      @author Anonymous (possibly A. Pedersen or G. Henkelman ), revision: Jean Claude C. Berthet
      @date Unknown, revision: 2010, University of Iceland
      */
//#include "LJBinary.h"
#include <cmath>
#include "common/PotentialsInterface.h"

class Morse : public PotentialsInterface {
public:
      Morse();
      Morse(double re, double De, double a, double cutoff);
      void cleanMemory(void);// required by PotentialsInterface
      void force(long N, const double *R, const long *, double *F, double *U, const double *box);
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
