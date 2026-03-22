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
#include "LammpsLoader.h"

#include <iostream>

namespace eonc {

LammpsLoader &LammpsLoader::instance() {
  static LammpsLoader loader;
  return loader;
}

LammpsLoader::LammpsLoader() {
#ifdef _WIN32
  const char *names[] = {"lammps.dll", "liblammps.dll", nullptr};
#elif defined(__APPLE__)
  const char *names[] = {"liblammps.dylib", "liblammps.0.dylib", nullptr};
#else
  const char *names[] = {"liblammps.so", "liblammps.so.0", nullptr};
#endif

  m_handle = dynlib::openFirst(names);
  if (!m_handle) {
    return; // Not found -- require_loaded() will throw if user requests LAMMPS
  }

  // Load required symbols
  open_no_mpi      = dynlib::loadSym<open_no_mpi_fn>(m_handle, "lammps_open_no_mpi");
  close            = dynlib::loadSym<close_fn>(m_handle, "lammps_close");
  command          = dynlib::loadSym<command_fn>(m_handle, "lammps_command");
  file             = dynlib::loadSym<file_fn>(m_handle, "lammps_file");
  scatter_atoms    = dynlib::loadSym<scatter_atoms_fn>(m_handle, "lammps_scatter_atoms");
  extract_variable = dynlib::loadSym<extract_var_fn>(m_handle, "lammps_extract_variable");

#ifdef EONMPI
  open_mpi = dynlib::loadSym<open_mpi_fn>(m_handle, "lammps_open");
#endif

  if (!open_no_mpi || !close || !command || !file ||
      !scatter_atoms || !extract_variable) {
    std::cerr << "[LAMMPS] Library loaded but missing required symbols\n";
    dynlib::close(m_handle);
    m_handle = {};
    return;
  }

  m_loaded = true;
}

LammpsLoader::~LammpsLoader() {
  dynlib::close(m_handle);
}

void LammpsLoader::require_loaded() const {
  if (!m_loaded) {
    throw std::runtime_error(
        "LAMMPS potential requested but liblammps not found.\n"
        "Install via: conda install -c conda-forge lammps\n"
        "Or ensure liblammps is in your library search path "
        "(LD_LIBRARY_PATH / DYLD_LIBRARY_PATH / PATH).");
  }
}

} // namespace eonc
