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
#include "client/Matter.h"
#include "client/io/IOHelpers.hpp"

#include <fstream>
#include <string>

namespace eonc::io {

class WriteBase {
public:
  virtual bool write(const Matter &, const std::string &, bool = false) = 0;
};

template <typename T> class Writer : public WriteBase {
public:
  // To be implemented by the child classes
  virtual bool writeImpl(const Matter &, std::ofstream &) = 0;

  bool write(const Matter &mat, const std::string &fname,
             bool append = false) override final {
    std::ofstream fout;
    ensureFileOpen(fout, fname, append);
    bool success = static_cast<T *>(this)->writeImpl(mat, fout);
    fout.close();
    return success;
  }
};

// Base class for reading
template <typename Derived> class ReadBase {
public:
  bool read(Matter &matter, const std::string &filename) {
    std::ifstream file;
    ensureFileOpen(file, filename);

    if (!file.is_open()) {
      return false;
    }

    bool success = static_cast<Derived *>(this)->readImpl(matter, file);
    file.close();
    return success;
  }
};

} // namespace eonc::io
