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
#include <cstdio>
#include <string>

namespace eonc {
class Matter;

namespace io {

// Reading
bool con2matter(Matter &m, std::string filename);
bool con2matter(Matter &m, FILE *file);
bool convel2matter(Matter &m, std::string filename);
bool convel2matter(Matter &m, FILE *file);

// Writing
bool matter2con(Matter &m, std::string filename, bool append = false);
bool matter2con(Matter &m, FILE *file);
bool matter2convel(Matter &m, std::string filename);
bool matter2convel(Matter &m, FILE *file);
void matter2xyz(Matter &m, std::string filename, bool append = false);
void writeTibble(Matter &m, std::string filename);

} // namespace io
} // namespace eonc
