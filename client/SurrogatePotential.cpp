#include "SurrogatePotential.h"

std::tuple<double, AtomMatrix, double>
SurrogatePotential::get_ef_var(const AtomMatrix pos, const Vector<int> atmnrs,
                               const Matrix3S box) {
  double energy{std::numeric_limits<double>::infinity()};
  long nAtoms{pos.rows()};
  AtomMatrix forces{MatrixType::Zero(nAtoms, 3)};
  double var{0};
  this->force(nAtoms, pos.data(), atmnrs.data(), forces.data(), &energy, &var,
              box.data());
  return std::make_tuple(energy, forces, var);
};
