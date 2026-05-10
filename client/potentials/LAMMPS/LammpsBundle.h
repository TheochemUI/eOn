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

/// LAMMPS portable run-input bundle (.eonlpb).
///
/// One self-contained blob holding everything liblammps needs to run:
/// in.lammps plus every pair_coeff data file (EAM, Tersoff, ZBL,
/// ReaxFF, ...), every read_data input, every dlopen-loaded
/// pair_style .so plugin, every KIM table, every shell helper. The
/// content is opaque to eOn -- the packer just walks a directory and
/// the client just untars it. eonclient extracts the bundle into a
/// private scratch dir under temp_directory_path() and issues
/// `shell cd <scratch>` to liblammps so every relative path inside
/// in.lammps resolves there. The eonclient CWD becomes irrelevant
/// for LAMMPS file lookups.
///
/// Format (all little-endian):
///   [0..7]   magic    : "EONLPB1\0"   (8 bytes)
///   [8..15]  m_len    : uint64_t      (manifest length in bytes)
///   [16..]   manifest : JSON UTF-8    (m_len bytes)
///   [...]    bodies   : concatenated file contents in manifest order
///
/// Manifest schema:
///   {"files":[{"name":"in.lammps","size":1234},
///             {"name":"Pd.eam.alloy","size":56789}]}
///
/// `name` is the relative path inside the bundle. After extraction
/// every relative reference inside in.lammps (pair_coeff <file>,
/// include <file>, read_data <file>, plugin load <file.so>, KIM tables,
/// ...) resolves against the scratch dir verbatim, because eOn does
/// `shell cd <scratch>` before `lammps_file("in.lammps")`.

#include <filesystem>
#include <string>
#include <vector>

namespace eonc {

struct LAMMPSBundleEntry {
  std::string name;
  std::uint64_t size;
};

class LAMMPSBundle {
public:
  /// Verify the magic and parse the manifest. Throws on a bad header
  /// or malformed manifest.
  static LAMMPSBundle open(const std::filesystem::path &bundle);

  /// Extract every entry to a fresh per-instance scratch dir under
  /// std::filesystem::temp_directory_path(). Returns the scratch dir
  /// path. The caller is responsible for std::filesystem::remove_all
  /// when done.
  std::filesystem::path extract() const;

  [[nodiscard]] const std::vector<LAMMPSBundleEntry> &entries() const noexcept {
    return m_entries;
  }

  /// Path that open() was given; the bundle is read again at extract()
  /// time so we don't pin the file in memory between calls.
  [[nodiscard]] const std::filesystem::path &source() const noexcept {
    return m_source;
  }

private:
  LAMMPSBundle() = default;

  std::filesystem::path m_source;
  std::vector<LAMMPSBundleEntry> m_entries;
  /// Byte offset where the first entry's body begins (= 16 + manifest len).
  std::uint64_t m_bodies_offset{0};
};

} // namespace eonc
