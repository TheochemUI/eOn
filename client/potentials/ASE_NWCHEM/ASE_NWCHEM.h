//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#pragma once
#define PYBIND11_DETAILED_ERROR_MESSAGES

#include "../../Potential.h"

#include <pybind11/eigen.h>
#include <pybind11/embed.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace pybind11::literals; // to bring in the `_a` literal

class ASENwchemPot : public Potential {

private:
  py::object calc;
  py::object ase;
  size_t counter;

public:
  ASENwchemPot(shared_ptr<Parameters> a_params);
  virtual ~ASENwchemPot() {
    SPDLOG_INFO("[ASENwchem] called potential {} times", counter);
  }

  // Functions
  void force(long nAtoms, const double *R, const int *atomicNrs,
             double *F, double *U, double *variance,
             const double *box) override;
};
