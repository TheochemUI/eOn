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
#include "OptimCreator.hpp"
#include "ConjugateGradients.h"

namespace eonc {
ConjugateGradients
mkOptim::operator()(const ObjectiveFunction &objf,
                    const ConjugateGradients::Params &params) {
  return ConjugateGradients(objf, params);
}
void mkOptim::operator()(const ObjectiveFunction &, const std::monostate &) {
  throw std::logic_error("Incorrect OptParams variant for OptBaseVisitor");
}
} // namespace eonc
