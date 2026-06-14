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

/// Runtime loader for Fortran potential shared libraries.
///
/// Uses dlopen (POSIX) or LoadLibrary (Windows) to load Fortran potential
/// libraries at runtime. The search order is:
///   1. Paths from the [Potential] potentials_path config key
///   2. EON_POTENTIALS_PATH environment variable (colon-separated)
///   3. Default system library search path (LD_LIBRARY_PATH, etc.)
///
/// Call add_config_paths() once at startup (before any potential constructor)
/// to inject config-file paths into the singleton.

#include "../DynLib.h"

#include <mutex>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace eonc {

class FortranPotLoader {
public:
  /// Thread-safe singleton accessor (Meyer's pattern).
  static FortranPotLoader &instance();

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

  /// Get the current search paths (for diagnostics).
  [[nodiscard]] const std::vector<std::string> &
  search_paths() const noexcept {
    return m_search_paths;
  }

  FortranPotLoader(const FortranPotLoader &) = delete;
  FortranPotLoader &operator=(const FortranPotLoader &) = delete;

private:
  FortranPotLoader();
  ~FortranPotLoader();

  dynlib::Handle open_lib(const char *lib_base);
  std::vector<std::string> lib_names(const char *lib_base) const;
  void append_paths(const std::string &path_str, char sep);

  std::vector<std::string> m_search_paths;
  std::mutex m_mutex;
  std::unordered_map<std::string, dynlib::Handle> m_handles;
};

} // namespace eonc
