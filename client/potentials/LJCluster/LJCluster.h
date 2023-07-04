//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef LENNARD_JONES_CLUSTER_H
#define LENNARD_JONES_CLUSTER_H

#include <iostream>
#include <math.h>
// #include "../../system_unit.h" // unit converters
#include "../../Potential.h"

/** Lennard Jones potential.*/
class LJCluster : public Potential {

private:
  //	Variables
  double u0;
  double cuttOffR;
  double psi;

  double cuttOffU;

public:
  // Functions
  // constructor
  LJCluster(std::shared_ptr<Parameters> params)
      : Potential(params), u0{1.0}, cuttOffR{15.0}, psi{1.0} {};
  // TODO: Put these in parameters
  // LJCluster(double r0Recieved, double u0Recieved, double psiRecieved);

  ~LJCluster();

  // Just to satisfy interface
  void cleanMemory();

  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box) override;
  void setParameters(double r0Recieved, double u0Recieved, double psiRecieved);
};
#endif
