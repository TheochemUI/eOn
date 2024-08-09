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
#include "client/Eigen.h"
#include <filesystem>

namespace eonc::io {
VectorType getUniqueValues(const VectorType &);
Vector<size_t> getUniqueCounts(const Vector<size_t> &);

void ensureFileOpen(std::ofstream &, const std::filesystem::path &, bool);
void ensureFileOpen(std::ifstream &, const std::filesystem::path &);
} // namespace eonc::io
