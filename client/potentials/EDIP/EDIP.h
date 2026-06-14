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
#include "../FortranPotLoader.h"

/// EDIP potential for Si (Fortran implementation, loaded at runtime).
class EDIP final : public Potential {
public:
  using force_fn = void (*)(const long int *N, const double *R, double *F,
                            double *U, const double *bx, const double *by,
                            const double *bz);

  explicit EDIP(const Parameters &params)
      : Potential(PotType::EDIP, params) {
    auto &loader = eonc::FortranPotLoader::instance();
    m_force = loader.load_sym<force_fn>("eon_edip", "edip_");
    if (!m_force) {
      loader.throw_not_found("eon_edip", "EDIP potential (Fortran)");
    }
  }
  ~EDIP() override = default;

  [[nodiscard]] bool isThreadSafe() const noexcept override { return false; }

  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box) override;

private:
  force_fn m_force{nullptr};
};
