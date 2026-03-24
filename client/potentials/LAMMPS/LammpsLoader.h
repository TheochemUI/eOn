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

/// Runtime loader for the LAMMPS C library.
///
/// Uses dlopen (POSIX) or LoadLibrary (Windows) to load liblammps at runtime
/// rather than requiring it at compile time. This allows a single eOn binary
/// to optionally use LAMMPS potentials if the library is installed.

#include "../../DynLib.h"

#include <stdexcept>
#include <string>

#ifdef EONMPI
#include <mpi.h>
#endif

namespace eonc {

// ---------------------------------------------------------------------------
// LammpsLoader: singleton that loads liblammps function pointers at runtime
// ---------------------------------------------------------------------------
class LammpsLoader {
public:
  // LAMMPS C API function pointer types (from library.h)
  using open_no_mpi_fn = void *(*)(int, char **, void **);
  using close_fn = void (*)(void *);
  using command_fn = char *(*)(void *, const char *);
  using file_fn = void (*)(void *, const char *);
  using scatter_atoms_fn = void (*)(void *, const char *, int, int, void *);
  using extract_var_fn = void *(*)(void *, const char *, const char *);
#ifdef EONMPI
  using open_mpi_fn = void *(*)(int, char **, MPI_Comm, void **);
#endif

  /// Thread-safe singleton accessor (Meyer's pattern).
  static LammpsLoader &instance();

  // Loaded function pointers (null if library not found)
  open_no_mpi_fn open_no_mpi{nullptr};
  close_fn close{nullptr};
  command_fn command{nullptr};
  file_fn file{nullptr};
  scatter_atoms_fn scatter_atoms{nullptr};
  extract_var_fn extract_variable{nullptr};
#ifdef EONMPI
  open_mpi_fn open_mpi{nullptr};
#endif

  /// True if liblammps was successfully loaded.
  [[nodiscard]] bool is_loaded() const noexcept { return m_loaded; }

  /// Throws std::runtime_error if liblammps is not available.
  void require_loaded() const;

  LammpsLoader(const LammpsLoader &) = delete;
  LammpsLoader &operator=(const LammpsLoader &) = delete;

private:
  LammpsLoader();
  ~LammpsLoader();

  bool m_loaded{false};
  dynlib::Handle m_handle{};

  /// Try to load a symbol; returns nullptr on failure.
  template <typename Fn> Fn load_sym(const char *name) {
    return reinterpret_cast<Fn>(dynlib::sym(m_handle, name));
  }
};

} // namespace eonc
