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
#include "LammpsBundle.h"

#include <nlohmann/json.hpp>

#include <cstdint>
#include <cstring>
#include <format>
#include <fstream>
#include <random>
#include <stdexcept>

#ifdef _WIN32
#include <process.h>
#define eonc_getpid _getpid
#else
#include <unistd.h>
#define eonc_getpid getpid
#endif

namespace eonc {

namespace {

constexpr char kMagic[8] = {'E', 'O', 'N', 'L', 'P', 'B', '1', '\0'};

std::uint64_t read_u64_le(const char *p) noexcept {
  std::uint64_t v = 0;
  std::memcpy(&v, p, sizeof(v));
  return v;
}

std::filesystem::path make_scratch_dir() {
  std::random_device rd;
  std::mt19937_64 rng{rd()};
  for (int attempt = 0; attempt < 16; ++attempt) {
    auto p = std::filesystem::temp_directory_path() /
             std::format("eonc_lammps_{}_{:016x}",
                         static_cast<long>(eonc_getpid()), rng());
    std::error_code ec;
    if (std::filesystem::create_directories(p, ec) && !ec) {
      return p;
    }
  }
  throw std::runtime_error(
      "LAMMPSBundle: failed to allocate scratch dir under " +
      std::filesystem::temp_directory_path().string());
}

} // namespace

LAMMPSBundle LAMMPSBundle::open(const std::filesystem::path &bundle) {
  if (!std::filesystem::exists(bundle)) {
    throw std::runtime_error("LAMMPSBundle: file not found: " +
                             bundle.string());
  }
  std::ifstream in(bundle, std::ios::binary);
  if (!in) {
    throw std::runtime_error("LAMMPSBundle: failed to open " + bundle.string());
  }

  char header[16]{};
  in.read(header, sizeof(header));
  if (in.gcount() != static_cast<std::streamsize>(sizeof(header))) {
    throw std::runtime_error("LAMMPSBundle: short read on header (" +
                             bundle.string() + ")");
  }
  if (std::memcmp(header, kMagic, sizeof(kMagic)) != 0) {
    throw std::runtime_error("LAMMPSBundle: bad magic (expected EONLPB1) in " +
                             bundle.string());
  }
  const auto manifest_len = read_u64_le(header + 8);
  if (manifest_len == 0 || manifest_len > (1ULL << 24)) {
    throw std::runtime_error(
        std::format("LAMMPSBundle: implausible manifest length {} in {}",
                    manifest_len, bundle.string()));
  }

  std::string manifest_text(manifest_len, '\0');
  in.read(manifest_text.data(), static_cast<std::streamsize>(manifest_len));
  if (in.gcount() != static_cast<std::streamsize>(manifest_len)) {
    throw std::runtime_error("LAMMPSBundle: short read on manifest (" +
                             bundle.string() + ")");
  }

  nlohmann::json j;
  try {
    j = nlohmann::json::parse(manifest_text);
  } catch (const std::exception &e) {
    throw std::runtime_error(
        std::string("LAMMPSBundle: manifest is not valid JSON: ") + e.what());
  }
  if (!j.contains("files") || !j["files"].is_array()) {
    throw std::runtime_error("LAMMPSBundle: manifest missing 'files' array (" +
                             bundle.string() + ")");
  }

  LAMMPSBundle out;
  out.m_source = bundle;
  out.m_bodies_offset = sizeof(header) + manifest_len;
  out.m_entries.reserve(j["files"].size());
  for (const auto &entry : j["files"]) {
    if (!entry.contains("name") || !entry.contains("size")) {
      throw std::runtime_error(
          "LAMMPSBundle: manifest entry missing name/size in " +
          bundle.string());
    }
    out.m_entries.push_back(
        {entry["name"].get<std::string>(), entry["size"].get<std::uint64_t>()});
  }
  return out;
}

std::filesystem::path LAMMPSBundle::extract() const {
  std::ifstream in(m_source, std::ios::binary);
  if (!in) {
    throw std::runtime_error("LAMMPSBundle: failed to reopen " +
                             m_source.string());
  }
  in.seekg(static_cast<std::streamoff>(m_bodies_offset));
  if (!in) {
    throw std::runtime_error("LAMMPSBundle: seek to bodies failed in " +
                             m_source.string());
  }

  auto scratch = make_scratch_dir();
  for (const auto &entry : m_entries) {
    // Reject path traversal: bundle entries must be plain names or
    // forward-slash relative paths that stay inside the scratch dir.
    std::filesystem::path rel(entry.name);
    if (rel.is_absolute() || entry.name.find("..") != std::string::npos) {
      std::error_code ec;
      std::filesystem::remove_all(scratch, ec);
      throw std::runtime_error("LAMMPSBundle: rejecting unsafe entry name '" +
                               entry.name + "' in " + m_source.string());
    }
    auto out_path = scratch / rel;
    std::filesystem::create_directories(out_path.parent_path());
    std::ofstream out(out_path, std::ios::binary);
    if (!out) {
      std::error_code ec;
      std::filesystem::remove_all(scratch, ec);
      throw std::runtime_error("LAMMPSBundle: cannot write " +
                               out_path.string());
    }

    constexpr std::size_t kBufSize = 64 * 1024;
    std::vector<char> buf(kBufSize);
    std::uint64_t remaining = entry.size;
    while (remaining > 0) {
      const auto chunk = std::min<std::uint64_t>(remaining, kBufSize);
      in.read(buf.data(), static_cast<std::streamsize>(chunk));
      if (in.gcount() != static_cast<std::streamsize>(chunk)) {
        std::error_code ec;
        std::filesystem::remove_all(scratch, ec);
        throw std::runtime_error(
            std::format("LAMMPSBundle: short read for '{}' ({} of {} bytes) "
                        "in {}",
                        entry.name, static_cast<std::uint64_t>(in.gcount()),
                        entry.size, m_source.string()));
      }
      out.write(buf.data(), static_cast<std::streamsize>(chunk));
      remaining -= chunk;
    }
  }
  return scratch;
}

} // namespace eonc
