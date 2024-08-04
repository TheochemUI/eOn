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
#include <magic_enum/magic_enum_all.hpp>

#include "client/io/IOBase.hpp"

namespace eonc::io {

enum class WriteType {
  UNKNOWN = -1, // error case
  NONE = 0,
  CON,
  CONVEL,
  XYZ,
  TIBBLE
};

std::unique_ptr<WriteBase> mkWriter(const toml::table &config);
} // namespace eonc::io
