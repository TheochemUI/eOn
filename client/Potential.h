#ifndef POTENTIAL_H
#define POTENTIAL_H

#include "Eigen.h"
#include "Parameters.h"
#include <algorithm>
#include <memory>

class Potential {
private:
    // should be const
  PotType ptype;
  Parameters *params;

public:
  Potential(Parameters *parameters) : ptype{parameters->potential} {}
  ~Potential() = default;

  static int fcalls;
  static int fcallsTotal;
  static int wu_fcallsTotal;
  static double totalUserTime;

  void virtual force(long nAtoms, const double *positions, const int *atomicNrs,
                     double *forces, double *energy, const double *box) = 0;
  PotType getType() { return this->ptype; };
};

namespace helper_functions {
Potential *makePotential(Parameters *params);
  // Does not take into account the fixed / free atoms
  // Free function for bindings
std::pair<double, AtomMatrix> efPot(Potential* pot, const AtomMatrix pos, const VectorXi atmnrs, const Matrix3d box);
}

// int Potential::fcalls = 0;
// int Potential::fcallsTotal = 0;
// int Potential::wu_fcallsTotal = 0;
// double Potential::totalUserTime = 0.0;
#endif
