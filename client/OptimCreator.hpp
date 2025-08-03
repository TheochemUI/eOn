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
#include "client/ConjugateGradients.h"
#include "client/ObjectiveFunction.h"
#include <variant>
namespace eonc {
using OptParams =
    std::variant<std::monostate, ConjugateGradients::Params /*, OtherParam */>;

struct mkOptim {
  ConjugateGradients operator()(const ObjectiveFunction &,
                                const ConjugateGradients::Params &);
  void operator()(const ObjectiveFunction &, const std::monostate &);
};

} // namespace eonc
