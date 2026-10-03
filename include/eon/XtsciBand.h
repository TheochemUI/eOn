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

#include "Eigen.h"
#include "Parameters.h"

#include <vector>

struct RgsaddleBand;

namespace eonc {

class NudgedElasticBand;

/// One rgsaddle band step. Matter owns positions; the session owns the
/// projected-force assembly and the inner stepper.
class XtsciBand {
public:
  XtsciBand(NudgedElasticBand &neb, const Parameters &params);
  ~XtsciBand();

  XtsciBand(const XtsciBand &) = delete;
  XtsciBand &operator=(const XtsciBand &) = delete;

  /// Copy the current path into the session.
  void syncFromPath();
  /// Drop stepper history after a host reparameterization.
  void reset();
  /// Assemble the band force and take one stepper step.
  void step(double maxMove);

private:
  friend int xtsciBandSurface(void *user, void *request);
  // True when the session's band is the host path, bit for bit.
  bool sessionMatchesPath() const;

  NudgedElasticBand *m_neb;
  RgsaddleBand *m_band{nullptr};
  std::vector<double> m_springKs;
  std::vector<double> m_cell;
  std::vector<AtomMatrix> m_fixed;
  double m_maxMove{0.2};
};

} // namespace eonc
