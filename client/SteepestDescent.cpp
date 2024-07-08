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

// Based on the SteepestDescent minimizer written in ASE.

#include "SteepestDescent.h"
#include "HelperFunctions.h"
namespace eonc {
int SteepestDescent::step(double a_maxMove) {
  VectorType r = m_objf->getPositions();
  VectorType f = -m_objf->getGradient();

  VectorType dr;
  double alpha = m_params->optim.SDAlpha;
  if (m_params->optim.SDTwoPoint == true && iteration > 0) {
    VectorType dx = r - m_rPrev;
    VectorType dg = -f + m_fPrev;
    alpha = dx.dot(dx) / dx.dot(dg);
    if (alpha < 0) {
      alpha = m_params->optim.SDAlpha;
    }
    SPDLOG_LOGGER_DEBUG(m_log, "[SD] alpha: {:.4e}", alpha);
  }

  dr = alpha * f;
  dr = helper_functions::maxAtomMotionAppliedV(dr, a_maxMove);

  m_objf->setPositions(r + dr);

  m_rPrev = r;
  m_fPrev = f;

  iteration++;

  return m_objf->isConverged() ? 1 : 0;
}

int SteepestDescent::run(size_t a_maxIteration, double a_maxMove) {
  while (!m_objf->isConverged() && iteration < a_maxIteration) {
    step(a_maxMove);
  }
  return m_objf->isConverged() ? 1 : 0;
}

} // namespace eonc
