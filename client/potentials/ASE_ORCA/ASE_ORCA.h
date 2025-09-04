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
#define PYBIND11_DETAILED_ERROR_MESSAGES

#include "../../Potential.h"

#include <pybind11/eigen.h>
#include <pybind11/embed.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace pybind11::literals; // to bring in the `_a` literal

class ASEOrcaPot : public Potential {

private:
  py::object calc;
  py::object ase;
  size_t counter;

public:
  ASEOrcaPot(shared_ptr<Parameters> a_params);
  virtual ~ASEOrcaPot() {
    SPDLOG_INFO("[ASEOrca] called potential {} times", counter);
  }

  // Functions
  void force(long nAtoms, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box) override;
};
