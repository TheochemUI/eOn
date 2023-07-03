#pragma once

#include "Potential.h"

class SurrogatePotential : public Potential {

public:
  SurrogatePotential(std::shared_ptr<Parameters> paramsIn)
      : Potential(paramsIn) {}
  virtual ~SurrogatePotential() = default;
  std::tuple<double, AtomMatrix, Eigen::MatrixXd>
  get_ef_var(const AtomMatrix pos, const VectorXi atmnrs, const Matrix3d box);
};
