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

/// Tersoff potential for Si (Fortran implementation, loaded at runtime).
class Tersoff final : public Potential {
public:
  using force_fn = void (*)(const long int *N, const double *R, double *F,
                            double *U, const double *bx, const double *by,
                            const double *bz);

  explicit Tersoff(const Parameters &p)
      : Potential(p) {
    auto &loader = eonc::FortranPotLoader::instance();
    m_force = loader.load_sym<force_fn>("eon_tersoff", "tersoff_");
    if (!m_force) {
      loader.throw_not_found("eon_tersoff", "Tersoff potential (Fortran)");
    }
  }
  ~Tersoff() override = default;

  [[nodiscard]] bool isThreadSafe() const noexcept override { return false; }

  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box) override;

private:
  force_fn m_force{nullptr};
};
