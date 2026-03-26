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
#include <array>
#include <readcon-core.hpp>
#include <string>
#include <utility>

namespace eonc {
class Matter;

namespace io {

// Reading
bool con2matter(Matter &m, std::string filename);
bool con2matter(Matter &m, const readcon::ConFrame &frame);
bool convel2matter(Matter &m, std::string filename);

// Writing
bool matter2con(Matter &m, std::string filename, bool append = false);
bool matter2convel(Matter &m, std::string filename);
void matter2xyz(Matter &m, std::string filename, bool append = false);
void writeTibble(Matter &m, std::string filename);

// Helper
std::pair<std::array<double, 3>, std::array<double, 3>>
cell_to_lengths_angles(const Matter &m);

} // namespace io
} // namespace eonc
