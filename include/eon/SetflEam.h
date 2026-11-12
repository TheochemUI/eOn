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

#include "Potential.h"

#include <memory>
#include <string>

namespace eonc {

/// Tabulated eam/alloy setfl potential with the Cauchy stress.
/// The Mishin Cu file (Phys. Rev. B 63, 224106) is the copper case.
class SetflEam final : public Potential {
public:
  SetflEam(const Parameters &params, std::string path);

  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box) override;

  [[nodiscard]] bool computesStress() const noexcept override { return true; }

  [[nodiscard]] Matrix3d cauchyStress() const override { return stress_; }

  /// The stress is stored on the instance, so two threads do not share it.
  [[nodiscard]] bool isThreadSafe() const noexcept override { return false; }

  [[nodiscard]] bool needsPerImageInstance() const noexcept override {
    return true;
  }

private:
  struct Tables;
  std::shared_ptr<const Tables> tables_;
  Matrix3d stress_{Matrix3d::Zero()};
};

} // namespace eonc
