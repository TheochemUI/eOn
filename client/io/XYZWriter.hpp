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
#include "client/io/IOBase.hpp"

namespace eonc::io {
class XYZWriter final : public Writer<XYZWriter> {
public:
  bool writeImpl(const Matter &mat, std::ofstream &fout) override;
};

} // namespace eonc::io
