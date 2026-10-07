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

/// True when the loaded rkrdb_open belongs to the linked readcon-db.
/// A missing library returns false and does not change IoStatus.
bool readcon_db_mirror_ok();

/// Drop a cached load so a later call reads EON_READCON_DB_LIBRARY again.
void readcon_db_mirror_reset();

/// Version of the loaded library. Empty when the load failed.
const char *readcon_db_loaded_version();

} // namespace eonc::io
