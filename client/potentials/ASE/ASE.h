/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
**
** Copyright (c) 2010--present, eOn Development Team
** All rights reserved.
**
** Repo:
** https://github.com/TheochemUI/eOn
*/
#pragma once

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
