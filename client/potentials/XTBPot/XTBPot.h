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
#include "units.hpp"

#include <cstddef>

class XTBPot final : public Potential {
public:
  XTBPot(const Parameters &p);
  ~XTBPot() override;

  // Disable copy to prevent double-free of Fortran pointers
  XTBPot(const XTBPot &) = delete;
  XTBPot &operator=(const XTBPot &) = delete;

  void cleanMemory(void);

  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box) override;

  /// XTB Fortran library uses per-instance state (env/calc).
  /// Thread-safe with separate instances; not safe on same instance.
  [[nodiscard]] bool isThreadSafe() const noexcept override { return false; }

  /// XTB restart.f90 uses global Fortran file units that collide when
  /// multiple environments exist in one process. Run sequentially with
  /// a single instance until upstream fixes unit management.
  [[nodiscard]] bool needsPerImageInstance() const noexcept override {
    return false;
  }

private:
  enum class GFNMethod { GFNFF, GFN0xTB, GFN1xTB, GFN2xTB };

  // Opaque libxtb handles. xtb.h declares them as
  // `typedef struct _xtb_TX* xtb_TX;`. We carry void* in the loader
  // world; the .cpp casts at the call sites where libxtb's types
  // matter.
  void *env{nullptr};
  void *calc{nullptr};
  void *mol{nullptr};
  void *res{nullptr};

  GFNMethod xtb_paramset;
  double xtb_acc;
  double xtb_electronic_temperature;
  std::size_t xtb_max_iter;
  double total_charge{0.0};
  int uhf{0};

  std::size_t counter{0};
  bool initialized{false};
};
