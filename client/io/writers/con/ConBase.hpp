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
#include "client/matter/Matter.h"

namespace eonc::io {
class ConBaseWriter {
public:
  bool writeBase(const Matter &, std::ofstream &, bool);

protected:
  virtual void writeHeader(const Matter &, std::ofstream &);
  virtual void writeAtoms(const Matter &, std::ofstream &,
                          const std::vector<size_t> &, const VectorType &);
  virtual void writeVelocities(const Matter &, std::ofstream &,
                               const std::vector<size_t> &, const VectorType &);

private:
  // The returned vector will store the starting index of each component
  std::vector<size_t> calculateComponentStartIndices(const Matter &, size_t);
};

} // namespace eonc::io
