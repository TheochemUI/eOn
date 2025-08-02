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
#include "client/Log.hpp" // IWYU pragma: keep
#include "thirdparty/toml.hpp"
namespace eonc {
toml::table commandLine(std::shared_ptr<spdlog::logger>, int, char **);
}
