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
#include "eon/potentials/LAMMPS/LammpsLoader.h"

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>
#include <system_error>

namespace eonc {

namespace {

const char *const *lammps_lib_names() noexcept {
#ifdef _WIN32
  static const char *names[] = {"lammps.dll", "liblammps.dll", nullptr};
#elif defined(__APPLE__)
  static const char *names[] = {"liblammps.dylib", "liblammps.0.dylib",
                                nullptr};
#else
  static const char *names[] = {"liblammps.so", "liblammps.so.0", nullptr};
#endif
  return names;
}

bool path_exists(const std::filesystem::path &p) {
  std::error_code ec;
  return std::filesystem::exists(p, ec);
}

bool lib_file_visible(const char *const *names) {
  for (const char *const *n = names; *n; ++n) {
    if (path_exists(*n))
      return true;
  }
#ifdef _WIN32
  const char *path = std::getenv("PATH");
  const char sep = ';';
#else
  const char *path = std::getenv("LD_LIBRARY_PATH");
#ifdef __APPLE__
  if (!path || !*path)
    path = std::getenv("DYLD_LIBRARY_PATH");
#endif
  const char sep = ':';
#endif
  if (!path || !*path)
    return false;
  const std::string paths(path);
  std::string::size_type start = 0;
  while (start < paths.size()) {
    auto pos = paths.find(sep, start);
    if (pos == std::string::npos)
      pos = paths.size();
    if (pos > start) {
      const std::filesystem::path dir(paths.substr(start, pos - start));
      for (const char *const *n = names; *n; ++n) {
        if (path_exists(dir / *n))
          return true;
      }
    }
    start = pos + 1;
  }
  return false;
}

} // namespace

LammpsLoader &LammpsLoader::instance() {
  static LammpsLoader loader;
  return loader;
}

void LammpsLoader::ensure_loaded() {
  if (m_tried)
    return;
  m_tried = true;

  m_handle = dynlib::openFirst(lammps_lib_names());
  if (!m_handle) {
    return;
  }

  open_no_mpi = dynlib::loadSym<open_no_mpi_fn>(m_handle, "lammps_open_no_mpi");
  close = dynlib::loadSym<close_fn>(m_handle, "lammps_close");
  command = dynlib::loadSym<command_fn>(m_handle, "lammps_command");
  file = dynlib::loadSym<file_fn>(m_handle, "lammps_file");
  scatter_atoms =
      dynlib::loadSym<scatter_atoms_fn>(m_handle, "lammps_scatter_atoms");
  extract_variable =
      dynlib::loadSym<extract_var_fn>(m_handle, "lammps_extract_variable");

#ifdef EONMPI
  open_mpi = dynlib::loadSym<open_mpi_fn>(m_handle, "lammps_open");
#endif

  if (!open_no_mpi || !close || !command || !file || !scatter_atoms ||
      !extract_variable) {
    std::cerr << "[LAMMPS] Library loaded but missing required symbols\n";
    dynlib::close(m_handle);
    m_handle = {};
    return;
  }

  m_loaded = true;
}

LammpsLoader::~LammpsLoader() { dynlib::close(m_handle); }

bool LammpsLoader::available() const {
  return lib_file_visible(lammps_lib_names());
}

void LammpsLoader::require_loaded() {
  ensure_loaded();
  if (!m_loaded) {
    throw std::runtime_error(
        "LAMMPS potential requested but liblammps not found.\n"
        "Install via: conda install -c conda-forge lammps\n"
        "Or ensure liblammps is in your library search path "
        "(LD_LIBRARY_PATH / DYLD_LIBRARY_PATH / PATH).");
  }
}

} // namespace eonc
