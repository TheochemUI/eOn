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

#include "eon/Potential.h"
#include <pybind11/pybind11.h>

namespace py = pybind11;

class ASE : public eonc::Potential {

private:
  size_t counter{0};
  py::module_ py_module; // Member to store the Python module
  py::object calculator; // Member to store the ASE calculator object
  py::object _calculate; // Member to store the Python function to calculate
                         // forces and energy
  py::object batch_calculate_;
  bool has_batch_{false};

public:
  ASE(const eonc::Parameters &a_params);
  virtual ~ASE() {
    QUILL_LOG_INFO(eonc::log::get(), "[ASE] called potential {} times",
                   counter);
  }

  void force(long nAtoms, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box) override;
  [[nodiscard]] bool isThreadSafe() const noexcept override { return false; }
  /// ASE calculators are independent objects. Most production calculators
  /// (VASP, Gaussian, ORCA, etc.) spawn subprocesses that release the GIL.
  [[nodiscard]] bool needsPerImageInstance() const noexcept override {
    return true;
  }
  [[nodiscard]] bool supportsBatchEvaluation() const noexcept override {
    return has_batch_;
  }
  void forceBatch(long nSystems, long nAtoms, const double *const *positions,
                  const int *const *atomicNrs, double *const *forces,
                  double *energies, double *variances,
                  const double *const *boxes) override;
};
