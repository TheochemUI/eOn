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

class ASENwchemPot : public Potential {

private:
  py::object calc;
  py::object ase;
  size_t counter{0};

public:
  ASENwchemPot(const Parameters &a_params);
  virtual ~ASENwchemPot() {
    QUILL_LOG_INFO(eonc::log::get(), "[ASENwchem] called potential {} times",
                   counter);
  }

  // Functions
  void force(long nAtoms, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box) override;
  [[nodiscard]] bool isThreadSafe() const noexcept override { return false; }
  /// NWChem runs as external subprocess; separate instances are independent.
  [[nodiscard]] bool needsPerImageInstance() const noexcept override {
    return true;
  }
  /// ASE-NWChem molecular SCF does not support PBC (#188).
  [[nodiscard]] bool requiresIsolatedMoleculeLayout() const noexcept override {
    return true;
  }
};
