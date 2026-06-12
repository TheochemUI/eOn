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
#include "FortranPotLoader.h"

#include <cstdlib>
#include <filesystem>
#include <sstream>

namespace eonc {

FortranPotLoader &FortranPotLoader::instance() {
  static FortranPotLoader loader;
  return loader;
}

FortranPotLoader::FortranPotLoader() {
  // Seed from EON_POTENTIALS_PATH env var (colon-separated on POSIX).
  // Config-file paths are prepended later via add_config_paths().
  const char *env = std::getenv("EON_POTENTIALS_PATH");
  if (env && *env) {
#ifdef _WIN32
    append_paths(env, ';');
#else
    append_paths(env, ':');
#endif
  }
}

FortranPotLoader::~FortranPotLoader() {
  for (auto &[name, handle] : m_handles) {
    dynlib::close(handle);
  }
}

void FortranPotLoader::append_paths(const std::string &path_str, char sep) {
  std::string::size_type start = 0;
  while (start < path_str.size()) {
    auto pos = path_str.find(sep, start);
    if (pos == std::string::npos)
      pos = path_str.size();
    if (pos > start) {
      std::string p = path_str.substr(start, pos - start);
      bool dup = false;
      for (const auto &existing : m_search_paths) {
        if (existing == p) {
          dup = true;
          break;
        }
      }
      if (!dup)
        m_search_paths.push_back(std::move(p));
    }
    start = pos + 1;
  }
}

void FortranPotLoader::add_config_paths(const std::string &colon_paths) {
  if (colon_paths.empty())
    return;
  std::lock_guard<std::mutex> lock(m_mutex);
  // Config paths have higher priority than env-var paths: insert at front.
  std::vector<std::string> new_paths;
#ifdef _WIN32
  const char sep = ';';
#else
  const char sep = ':';
#endif
  std::string::size_type start = 0;
  while (start < colon_paths.size()) {
    auto pos = colon_paths.find(sep, start);
    if (pos == std::string::npos)
      pos = colon_paths.size();
    if (pos > start) {
      std::string p = colon_paths.substr(start, pos - start);
      bool dup = false;
      for (const auto &existing : m_search_paths) {
        if (existing == p) {
          dup = true;
          break;
        }
      }
      if (!dup)
        new_paths.push_back(std::move(p));
    }
    start = pos + 1;
  }
  m_search_paths.insert(m_search_paths.begin(), new_paths.begin(),
                        new_paths.end());
}

std::vector<std::string>
FortranPotLoader::lib_names(const char *lib_base) const {
  std::vector<std::string> names;
#ifdef _WIN32
  names.push_back(std::string(lib_base) + ".dll");
  names.push_back(std::string("lib") + lib_base + ".dll");
#elif defined(__APPLE__)
  names.push_back(std::string("lib") + lib_base + ".dylib");
#else
  names.push_back(std::string("lib") + lib_base + ".so");
#endif
  return names;
}

dynlib::Handle FortranPotLoader::open_lib(const char *lib_base) {
  std::lock_guard<std::mutex> lock(m_mutex);

  auto it = m_handles.find(lib_base);
  if (it != m_handles.end()) {
    return it->second;
  }

  auto names = lib_names(lib_base);
  dynlib::Handle h{};

  // Try each search path (config paths precede env-var paths in m_search_paths)
  for (const auto &dir : m_search_paths) {
    for (const auto &name : names) {
      std::string full = (std::filesystem::path(dir) / name).string();
      h = dynlib::open(full.c_str());
      if (h) {
        m_handles[lib_base] = h;
        return h;
      }
    }
  }

  // Fall back to system search (LD_LIBRARY_PATH, etc.)
  for (const auto &name : names) {
    h = dynlib::open(name.c_str());
    if (h) {
      m_handles[lib_base] = h;
      return h;
    }
  }

  m_handles[lib_base] = nullptr;
  return nullptr;
}

void FortranPotLoader::throw_not_found(const char *lib_base,
                                       const char *description) const {
  auto names = lib_names(lib_base);
  std::ostringstream oss;
  oss << description << " requested but library not found.\n"
      << "Looked for: ";
  for (size_t i = 0; i < names.size(); ++i) {
    if (i > 0)
      oss << ", ";
    oss << names[i];
  }
  oss << "\n";
  if (!m_search_paths.empty()) {
    oss << "Search paths (config potentials_path + EON_POTENTIALS_PATH):\n";
    for (const auto &p : m_search_paths) {
      oss << "  " << p << "\n";
    }
  }
  oss << "Set [Potential] potentials_path in config.ini, "
         "or set EON_POTENTIALS_PATH,\n"
         "or ensure the libraries are in LD_LIBRARY_PATH.";
  throw std::runtime_error(oss.str());
}

} // namespace eonc
