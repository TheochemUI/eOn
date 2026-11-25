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

#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>

namespace eonc::io {

/// FNV-1a 64 of the canonical path bytes. No NUL and no file text.
/// Appended frames extend this trajectory. A zero digest becomes 1.
/// Shared with eon/concorpus.py traj_id.
[[nodiscard]] std::uint64_t corpus_path_id(std::string_view canonical_path);

/// FNV-1a 64 of the canonical path, a NUL, and the con text. A whole-file
/// blob is keyed by content because the corpus has no replace: a rewrite
/// stores a new blob. Shared with eon/concorpus.py blob_id.
[[nodiscard]] std::uint64_t corpus_blob_id(std::string_view canonical_path,
                                           std::string_view text);

/// Copy a con file into the run readcon-db corpus when libreadcon_db loads.
/// An overwrite (w='w') inserts immediately. Failure to load the library or
/// to insert the blob does not change IoStatus.
void mirror_con_corpus(const std::string &path);

/// True when the loaded rkrdb_open belongs to the linked readcon-db.
/// A missing library returns false and does not change IoStatus.
bool readcon_db_mirror_ok();

/// Drop a cached load so a later call reads EON_READCON_DB_LIBRARY again.
void readcon_db_mirror_reset();

/// Version of the loaded library. Empty when the load failed.
const char *readcon_db_loaded_version();

/// Buffer one serialized frame from an append. At 64 frames the batch is
/// extended once. A missing extend symbol soft-skips.
void note_corpus_append(const std::string &path, const std::string &frame_text);

/// True when the loaded readcon-db can extend a trajectory, so a
/// buffered append batch has somewhere to go.
[[nodiscard]] bool corpus_extend_ok();

/// Extend every buffered batch. Called on corpus-cache teardown.
void flush_corpus_appends();

/// Frames waiting for the next extend of this path. Zero when the mirror
/// is compiled out.
[[nodiscard]] std::size_t corpus_buffered_frames(const std::string &path);

} // namespace eonc::io
