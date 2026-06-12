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

/// Aluminum EAM potential (Fortran implementation, loaded at runtime).
class Aluminum final : public Potential {
public:
  using force_fn = void (*)(const long int *N, const double *R, double *F,
                            double *U, const double *bx, const double *by,
                            const double *bz);
  using potinit_fn = void (*)();

  explicit Aluminum(const Parameters &params)
      : Potential(PotType::EAM_AL, params) {
    auto &loader = eonc::FortranPotLoader::instance();
    m_force = loader.load_sym<force_fn>("eon_aluminum", "force_");
    m_potinit = loader.load_sym<potinit_fn>("eon_aluminum", "potinit_");
    if (!m_force || !m_potinit) {
      loader.throw_not_found("eon_aluminum",
                             "Aluminum EAM potential (Fortran)");
    }
    m_potinit();
  }
  ~Aluminum() override = default;

  [[nodiscard]] bool isThreadSafe() const noexcept override;

  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box) override;

private:
  force_fn m_force{nullptr};
  potinit_fn m_potinit{nullptr};
};
