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
#include "client/io/writers/con/ConBase.hpp"

namespace eonc::io {
class ConWriter final : public Writer<ConWriter>, public ConBaseWriter {
public:
  bool writeImpl(const Matter &, std::ofstream &) override;
};

} // namespace eonc::io
