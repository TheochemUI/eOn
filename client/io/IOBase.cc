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
#include "client/io/IOBase.hpp"
#include <filesystem>

namespace eonc::io {

bool IOBase::fileExists(const std::string &filename) {
  return std::filesystem::exists(filename);
}

void IOBase::ensureFileOpen(std::ofstream &file, const std::string &filename,
                            bool append) {
  if (append) {
    file.open(filename, std::ios::out | std::ios::app);
  } else {
    file.open(filename, std::ios::out | std::ios::trunc);
  }
}

void IOBase::ensureFileOpen(std::ifstream &file, const std::string &filename) {
  file.open(filename, std::ios::in);
}

} // namespace eonc::io
