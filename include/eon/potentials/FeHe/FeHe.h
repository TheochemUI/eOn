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
#include "eon/potentials/FortranPotLoader.h"

/// Fe-He interatomic potential (Fortran implementation, loaded at runtime).
class FeHe final : public Potential {
public:
  using feforce_fn = void (*)(const long int *N, const double *RX,
                              const double *RY, const double *RZ,
                              const int *ISPEC, double *FX, double *FY,
                              double *FZ, double *U, const double *bx,
                              const double *by, const double *bz);

  explicit FeHe(const Parameters &params)
      : Potential(params) {
    auto &loader = eonc::FortranPotLoader::instance();
    m_feforce = loader.load_sym<feforce_fn>("eon_fehe", "feforce_");
    if (!m_feforce) {
      loader.throw_not_found("eon_fehe", "Fe-He potential (Fortran)");
    }
  }
  ~FeHe() override = default;

  [[nodiscard]] bool isThreadSafe() const noexcept override { return false; }
  /// Fortran backend keeps state in module/COMMON storage; a shared
  /// instance across threads is never safe.
  [[nodiscard]] bool isSharedInstanceThreadSafe() const noexcept override {
    return false;
  }

  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box) override;

private:
  feforce_fn m_feforce{nullptr};
};
