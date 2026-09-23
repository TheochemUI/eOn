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

#include <string>

namespace eonc::io {

/// Copy a con file into the run readcon-db corpus when libreadcon_db loads.
/// Failure to load the library or to insert the blob does not change IoStatus.
void mirror_con_corpus(const std::string &path);

} // namespace eonc::io
