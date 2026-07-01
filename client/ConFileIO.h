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
#include <array>
#include <cstdint>
#include <memory>
#include <optional>
#include <readcon-core.hpp>
#include <string>
#include <utility>
#include <vector>

namespace eonc {
class Matter;

namespace io {

/// Structured I/O result for the client surface (nanobind-friendly).
/// Prefer comparing to IoStatus::Ok rather than treating as bool.
enum class IoStatus : std::uint8_t {
  Ok = 0,
  ReadError = 1,
  WriteError = 2,
  AppendError = 3,
  OpenError = 4,
  InvalidArgument = 5,
};

[[nodiscard]] constexpr bool io_ok(IoStatus s) noexcept {
  return s == IoStatus::Ok;
}

/// Human-readable label for logging / bindings.
[[nodiscard]] constexpr const char *io_status_name(IoStatus s) noexcept {
  switch (s) {
  case IoStatus::Ok:
    return "ok";
  case IoStatus::ReadError:
    return "read_error";
  case IoStatus::WriteError:
    return "write_error";
  case IoStatus::AppendError:
    return "append_error";
  case IoStatus::OpenError:
    return "open_error";
  case IoStatus::InvalidArgument:
    return "invalid_argument";
  }
  return "unknown";
}

struct ConMetadataValue {
  std::string key;
  double value;
};

struct ConMetadataText {
  std::string key;
  std::string value;
};

struct ConFrameMetadata {
  std::optional<uint64_t> frame_index;
  std::optional<double> energy;
  std::optional<double> time;
  std::optional<double> timestep;
  std::optional<uint64_t> neb_bead;
  std::optional<uint64_t> neb_band;
  std::optional<std::string> potential_type;
  std::vector<ConMetadataValue> scalars;
  std::vector<ConMetadataText> strings;
  std::optional<std::string> raw_json;
};

/// Extract known frame-level fields from a parsed readcon frame.
ConFrameMetadata metadata_from_frame(const readcon::ConFrame &frame);

// Reading
[[nodiscard]] IoStatus con2matter(Matter &m, std::string filename);
[[nodiscard]] IoStatus con2matter(Matter &m, const readcon::ConFrame &frame,
                                  ConFrameMetadata *out_metadata = nullptr);
[[nodiscard]] IoStatus convel2matter(Matter &m, std::string filename);

// Writing
[[nodiscard]] IoStatus matter2con(Matter &m, std::string filename,
                                  bool append = false,
                                  const ConFrameMetadata *metadata = nullptr);
[[nodiscard]] IoStatus matter2convel(Matter &m, std::string filename);
[[nodiscard]] IoStatus matter2xyz(Matter &m, std::string filename,
                                  bool append = false);
[[nodiscard]] IoStatus writeTibble(Matter &m, std::string filename);

/**
 * Write a full NEB path as one multi-frame .con using ConFrameBuilder::clone().
 *
 * Seeds identity (symbols/fixed/mass/id/cell headers) from path[0] only, then
 * for each image clones that template and bulk-updates positions/forces +
 * metadata. Images must share atom count and topology with path[0] (standard
 * NEB band invariant); heterogeneous multi-frame movies should not use this.
 * Avoids re-reading the output file per image (legacy append path).
 *
 * @param path length must be numImages+2 (endpoints included)
 * @param metadata_per_image length must equal path.size()
 */
[[nodiscard]] IoStatus
writeNebPath(std::string filename,
             const std::vector<std::shared_ptr<Matter>> &path,
             const std::vector<ConFrameMetadata> &metadata_per_image);

// Helper
std::pair<std::array<double, 3>, std::array<double, 3>>
cell_to_lengths_angles(const Matter &m);

} // namespace io
} // namespace eonc
