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

#include "client/MinimizationJob.hpp"
#include "client/matter/StructComparer.hpp"

namespace eonc::job {
void from_toml(mat::StructComparer::Params &,
               const toml::node_view<const toml::node> &);
void from_toml(MinimizationJob::Params &,
               const toml::node_view<const toml::node> &);
} // namespace eonc::job
