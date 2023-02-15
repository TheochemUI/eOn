#ifndef POTENTIAL_H
#define POTENTIAL_H

#include "Eigen.h"
#include "Parameters.h"
#include <algorithm>
#include <memory>

class Potential {
private:
  PotType ptype;
  Parameters *params;

public:
  Potential(Parameters *parameters) : ptype{parameters->potential} {}
  ~Potential() = default;

  static int fcalls;
  static int fcallsTotal;
  static int wu_fcallsTotal;
  static double totalUserTime;

  void virtual initialize() = 0;
  void virtual force(long nAtoms, const double *positions, const int *atomicNrs,
                     double *forces, double *energy, const double *box) = 0;
  std::pair<double, AtomMatrix> get_ef(long nAtoms, const double *positions, const int *atomicNrs,
                     const double *box);
};

namespace helper_functions {
Potential *makePotential(Parameters *params);
}

// int Potential::fcalls = 0;
// int Potential::fcallsTotal = 0;
// int Potential::wu_fcallsTotal = 0;
// double Potential::totalUserTime = 0.0;
#endif
