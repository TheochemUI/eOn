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
#include "client/Optimizer.h"

namespace eonc::opt {
void from_toml(OptimBase::Params &, const toml::node_view<const toml::node> &);
void from_toml(ConjugateGradients::Params &,
               const toml::node_view<const toml::node> &);
} // namespace eonc::opt
