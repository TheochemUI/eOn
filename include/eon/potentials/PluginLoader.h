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

/// Runtime loader for engine plugin shared libraries.
///
/// Engines that ship as their own shared objects (the rgpot metatomic and
/// xtb backends) are found by base name across a search path rather than
/// by an absolute path, so a build, a wheel, and a conda prefix can each
/// place them differently.
///
/// Uses dlopen (POSIX) or LoadLibrary (Windows) to load Fortran potential
/// libraries at runtime. The search order is:
///   1. Paths from the [Potential] potentials_path config key
///   2. EON_POTENTIALS_PATH environment variable (';' on Windows, ':'
///   elsewhere)
///   3. Default system library search path (LD_LIBRARY_PATH, etc.)
///
/// Call add_config_paths() once at startup (before any potential constructor)
/// to inject config-file paths into the singleton.

#include "eon/DynLib.h"

#include <mutex>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace eonc {

class PluginLoader {
public:
  /// Thread-safe singleton accessor (Meyer's pattern).
  static PluginLoader &instance();

  /// Inject search paths from the eOn config file.
  /// @param colon_paths  colon-separated directory list (may be empty)
  /// These paths are inserted before the EON_POTENTIALS_PATH paths.
  void add_config_paths(const std::string &colon_paths);

  /// Load a symbol from a named potential library.
  /// @param lib_base  Base name without platform prefix/suffix (e.g. "eon_sw")
  /// @param sym_name  Symbol to load (e.g. "sw_")
  /// @returns Function pointer, or nullptr if not found.
  template <typename Fn>
  Fn load_sym(const char *lib_base, const char *sym_name) {
    dynlib::Handle h = open_lib(lib_base);
    if (!h)
      return nullptr;
    return dynlib::loadSym<Fn>(h, sym_name);
  }

  /// Throw a descriptive error for a missing potential library.
  [[noreturn]] void throw_not_found(const char *lib_base,
                                    const char *description) const;

  /// Filesystem probe: a matching library file is on the configured
  /// search path. Does not dlopen, so plugin static initializers
  /// (banners) do not run.
  [[nodiscard]] bool lib_present(const char *lib_base) const;

  /// Get the current search paths (for diagnostics).
  [[nodiscard]] const std::vector<std::string> &search_paths() const noexcept {
    return m_search_paths;
  }

  PluginLoader(const PluginLoader &) = delete;
  PluginLoader &operator=(const PluginLoader &) = delete;

private:
  PluginLoader();
  ~PluginLoader();

  dynlib::Handle open_lib(const char *lib_base);
  std::vector<std::string> lib_names(const char *lib_base) const;
  void append_paths(const std::string &path_str, char sep);

  std::vector<std::string> m_search_paths;
  mutable std::mutex m_mutex;
  std::unordered_map<std::string, dynlib::Handle> m_handles;
  /// Last dynlib::error() from a failed open (for throw_not_found).
  std::string m_last_load_error;
  bool m_last_load_saw_file = false;
};

} // namespace eonc
