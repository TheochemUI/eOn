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
#include <functional>
#include <string>
#include <string_view>
#include <vector>

namespace eonc::io {

/// FNV-1a 64 of the canonical path, a NUL, and the con text. A whole-file
/// blob is keyed by content because the corpus has no replace: a rewrite
/// stores a new blob. Shared with eon/concorpus.py blob_id.
[[nodiscard]] std::uint64_t corpus_blob_id(std::string_view canonical_path,
                                           std::string_view text);

/// Trajectory id of one generation of a path whose frame 0 is first_text.
/// Generation 0 is corpus_blob_id, so the blob an overwrite stores is the
/// trajectory its appends extend. Later generations hash the path, a NUL,
/// the text, a NUL and the decimal generation. A zero digest becomes 1.
/// Shared with eon/concorpus.py generation_id.
[[nodiscard]] std::uint64_t
corpus_generation_id(std::string_view canonical_path,
                     std::string_view first_text, std::uint64_t generation = 0);

/// Start the generation an overwrite of path with the one-frame text
/// begins. Frames buffered for the previous generation go to it first.
/// A generation that holds only this frame is reused.
void start_corpus_generation(const std::string &path, const std::string &text);

/// Drop the path's generation after a write that starts none, such as a
/// multi-frame snapshot. Buffered frames go to the old generation first.
void forget_corpus_generation(const std::string &path);

/// Per-frame con texts of a file, read only when a generation is resumed.
using CorpusDiskFrames = std::function<std::vector<std::string>()>;

/// Copy a con file into the run readcon-db corpus when libreadcon_db loads.
/// Reads and multi-frame snapshots insert immediately.
/// Failure to load the library or to insert the blob does not change IoStatus.
void mirror_con_corpus(const std::string &path);

/// True when the loaded rkrdb_open belongs to the linked readcon-db.
/// A missing library returns false and does not change IoStatus.
bool readcon_db_mirror_ok();

/// Drop a cached load so a later call reads EON_READCON_DB_LIBRARY again.
void readcon_db_mirror_reset();

/// Version of the loaded library. Empty when the load failed.
const char *readcon_db_loaded_version();

/// Buffer one serialized frame from an append for the path's generation.
/// At 64 frames the batch is extended once. A missing extend symbol
/// soft-skips. A path with no generation in this process resumes the
/// newest one for the file's first frame and brings it up to disk_frames,
/// which already hold this frame; stored frames that disagree with the
/// file start a new generation.
void note_corpus_append(const std::string &path, const std::string &frame_text,
                        const CorpusDiskFrames &disk_frames);

/// True when the loaded readcon-db can extend a trajectory, so a
/// buffered append batch has somewhere to go.
[[nodiscard]] bool corpus_extend_ok();

/// Extend every buffered batch. Called on corpus-cache teardown.
void flush_corpus_appends();

/// Frames waiting for the next extend of this path. Zero when the mirror
/// is compiled out.
[[nodiscard]] std::size_t corpus_buffered_frames(const std::string &path);

/// Frame texts of the path's current generation after a flush. Empty when
/// this process holds no generation for the path.
[[nodiscard]] std::vector<std::string>
corpus_generation_frame_texts(const std::string &path);

} // namespace eonc::io
