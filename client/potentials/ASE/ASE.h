//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef ASE_POTENTIAL
#define ASE_POTENTIAL

#include "../../Potential.h"
#include <pybind11/embed.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

class ASE : public Potential {

private:
  size_t counter;
  // XXX(rg): Remove this, guards are done in the client
  py::scoped_interpreter guard; // Member to manage the Python Interpreter
  py::module_ py_module;        // Member to store the Python module
  py::object calculator;        // Member to store the ASE calculator object
  py::object _calculate; // Member to store the Python function to calculate
                         // forces and energy

public:
  ASE(shared_ptr<Parameters> a_params);
  virtual ~ASE() { SPDLOG_INFO("[ASE] called potential {} times", counter); }

  void force(long nAtoms, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box) override;
};
#endif
