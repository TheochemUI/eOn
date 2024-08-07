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
#include "client/matter/Matter.h"
#include <fstream>
#include <optional>

namespace eonc::mat {

std::optional<std::string> readLine(std::ifstream &);
std::array<double, 3> readArray(const std::string &);

class ConFileParser {
private:
  size_t Ncomponent;
  std::array<std::string, 4> headers;
  // Vector of indices tracking first occurrences
  std::vector<size_t> first;

public:
  ConFileParser() = default;
  void parse(Matter &, const std::string &);

private:
  void con_headers(std::ifstream &, Matter &);
  void read_cartesian_fix_index(std::ifstream &, Matter &, AtomMatrix &,
                                const bool set_atmnr_);
  void finalize_matter(Matter &);
};

void from_con(Matter &, const std::string &);
void from_convel(Matter &, const std::string &);

void conCell(Matter &, const std::array<double, 3> &,
             const std::array<double, 3> &);

} // namespace eonc::mat
