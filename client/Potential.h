#ifndef POTENTIAL_H
#define POTENTIAL_H

#include "Eigen.h"
#include "Parameters.h"
#include <algorithm>
#include <limits>
#include <memory>

class Potential {
private:
  // should be const
  PotType ptype;

protected:
  std::shared_ptr<Parameters> params;

public:
  Potential(std::shared_ptr<Parameters> paramsIn)
      : ptype{paramsIn->potential} {}
  virtual ~Potential() = default;

  static int fcalls;
  static int fcallsTotal;
  static int wu_fcallsTotal;
  static double totalUserTime;

  // Does not take into account the fixed / free atoms
  void virtual force(long nAtoms, const double *positions, const int *atomicNrs,
                     double *forces, double *energy, const double *box) = 0;
  // Only implemented for surrogates
  // TODO: Maybe a subclass of Potential (SurrogatePotential) and use that
    std::tuple<double, AtomMatrix, Eigen::MatrixXd>
    get_ef_var(const AtomMatrix pos, const VectorXi atmnrs, const Matrix3d box);
  std::pair<double, AtomMatrix>
  get_ef(const AtomMatrix pos, const VectorXi atmnrs, const Matrix3d box) {
    double energy{std::numeric_limits<double>::infinity()};
    long nAtoms{pos.rows()};
    AtomMatrix forces{Eigen::MatrixXd::Zero(nAtoms, 3)};
    this->force(nAtoms, pos.data(), atmnrs.data(), forces.data(), &energy,
                box.data());
    return std::make_pair(energy, forces);
  };
  PotType getType() { return this->ptype; };
};

namespace helper_functions {
std::shared_ptr<Potential> makePotential(std::shared_ptr<Parameters> params);
std::shared_ptr<Potential> makePotential(PotType ptype,
                                         std::shared_ptr<Parameters> params);
} // namespace helper_functions

// int Potential::fcalls = 0;
// int Potential::fcallsTotal = 0;
// int Potential::wu_fcallsTotal = 0;
// double Potential::totalUserTime = 0.0;
#endif
