#pragma once
#include "Eigen.h"

#include "Matter.h"
#include "Parameters.h"

class Hessian {
public:
  Hessian(Parameters *params, Matter *matter);
  ~Hessian();

  MatrixType getHessian(Matter *matterIn, Vector<int> atomsIn);
  VectorType getFreqs(Matter *matterIn, Vector<int> atomsIn);
  //    VectorType getModes(Matter *matterIn, Vector<int> atomsIn);
  VectorType removeZeroFreqs(VectorType freqs);

private:
  Matter *matter;
  Parameters *parameters;

  MatrixType hessian;
  //    VectorType modes;
  VectorType freqs;

  Vector<int> atoms;
  bool calculate(void);
  std::shared_ptr<spdlog::logger> log;
};
