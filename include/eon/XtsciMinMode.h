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

#include "LowestEigenmode.h"

namespace eonc {

/// Lowest-mode rotation through a rgsaddle minimum-mode session.
/// The session is created with a zero translation cap so the climb
/// stays in MinModeSaddleSearch.
class XtsciMinMode : public LowestEigenmode {
public:
  XtsciMinMode(std::shared_ptr<Matter> matter, const Parameters &params,
               std::shared_ptr<Potential> pot);
  ~XtsciMinMode() override = default;

  void compute(std::shared_ptr<Matter> matter,
               AtomMatrix initialDirection) override;
  double getEigenvalue() override;
  AtomMatrix getEigenvector() override;

private:
  double m_eigenvalue{0.0};
  AtomMatrix m_eigenvector;
  AtomMatrix m_fixed;
};

} // namespace eonc
