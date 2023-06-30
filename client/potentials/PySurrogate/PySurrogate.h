//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef PYSURROGATE_INTERFACE
#define PYSURROGATE_INTERFACE

#define PYBIND11_DETAILED_ERROR_MESSAGES

#include "../../Potential.h"

#include <pybind11/eigen.h>
#include <pybind11/embed.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace pybind11::literals; // to bring in the `_a` literal

// TODO: Use a better name, say, CatLearnPot
class PySurrogate : public Potential {

private:
  py::object hpfit;
  py::object kernel;

public:
  PySurrogate(shared_ptr<Parameters> p);

  // Functions
  void train_optimize(Eigen::MatrixXd features, Eigen::MatrixXd targets);
  // To satisfy interface
  void cleanMemory(void);
  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box) override;

  // Variables [public]
  py::object gpmod;
  Eigen::MatrixXd variance; // XXX: This is a hacky way to populate and use this variable
};
#endif
