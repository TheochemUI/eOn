//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef CATLEARNPOT_INTERFACE
#define CATLEARNPOT_INTERFACE

#define PYBIND11_DETAILED_ERROR_MESSAGES

#include "../../SurrogatePotential.h"

#include <pybind11/eigen.h>
#include <pybind11/embed.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace pybind11::literals; // to bring in the `_a` literal

class CatLearnPot final : public SurrogatePotential {

public:
  CatLearnPot(shared_ptr<Parameters> a_params);
  py::object m_gpmod;

  // Functions
  void train_optimize(Eigen::MatrixXd a_features,
                      Eigen::MatrixXd a_targets) override;
  void force(long nAtoms, const double *positions, const int *atomicNrs,
             double *forces, double *energy, double *variance,
             const double *box) override;
  // Variables [public]
  py::object m_gpmod;
  Eigen::MatrixXd
      variance; // XXX: This is a hacky way to populate and use this variable
};
#endif
