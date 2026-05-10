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
#include <optional>
#include <readcon-core.hpp>
#include <string>
#include <utility>
#include <vector>

namespace eonc {
class Matter;

namespace io {

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
  std::vector<ConMetadataValue> scalars;
  std::vector<ConMetadataText> strings;
  std::optional<std::string> raw_json;
};

// Reading
bool con2matter(Matter &m, std::string filename);
bool con2matter(Matter &m, const readcon::ConFrame &frame);
bool convel2matter(Matter &m, std::string filename);

// Writing
bool matter2con(Matter &m, std::string filename, bool append = false,
                const ConFrameMetadata *metadata = nullptr);
bool matter2convel(Matter &m, std::string filename);
void matter2xyz(Matter &m, std::string filename, bool append = false);
void writeTibble(Matter &m, const std::string &filename);

// Helper
std::pair<std::array<double, 3>, std::array<double, 3>>
cell_to_lengths_angles(const Matter &m);

} // namespace io
} // namespace eonc
