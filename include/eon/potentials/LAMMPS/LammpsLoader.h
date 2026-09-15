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

#include "eon/DynLib.h"

#include <stdexcept>
#include <string>

#ifdef EONMPI
#include <mpi.h>
#endif

namespace eonc {

// ---------------------------------------------------------------------------
// ILammpsLoader: injectable ABI. Production uses LammpsLoader::instance().
// ---------------------------------------------------------------------------
class ILammpsLoader {
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

  open_no_mpi_fn open_no_mpi{nullptr};
  close_fn close{nullptr};
  command_fn command{nullptr};
  file_fn file{nullptr};
  scatter_atoms_fn scatter_atoms{nullptr};
  extract_var_fn extract_variable{nullptr};
#ifdef EONMPI
  open_mpi_fn open_mpi{nullptr};
#endif

  virtual ~ILammpsLoader() = default;
  virtual void require_loaded() = 0;
  [[nodiscard]] virtual bool is_loaded() const noexcept = 0;
  [[nodiscard]] virtual bool available() const = 0;
  [[nodiscard]] virtual const std::string &last_error() const noexcept = 0;

  ILammpsLoader(const ILammpsLoader &) = delete;
  ILammpsLoader &operator=(const ILammpsLoader &) = delete;

protected:
  ILammpsLoader() = default;
};

// ---------------------------------------------------------------------------
// LammpsLoader: process-default loader (Meyer's singleton). Does not dlopen
// until require_loaded().
// ---------------------------------------------------------------------------
class LammpsLoader : public ILammpsLoader {
public:
  /// Thread-safe singleton accessor (Meyer's pattern). Does not dlopen.
  static LammpsLoader &instance();

  [[nodiscard]] bool is_loaded() const noexcept override { return m_loaded; }

  /// Filesystem probe: a candidate liblammps file is visible without
  /// dlopen, so LAMMPS static initializers (the banner) do not run.
  [[nodiscard]] bool available() const override;

  /// Why the last ensure_loaded() failed. Empty if unused or loaded.
  [[nodiscard]] const std::string &last_error() const noexcept override {
    return m_last_error;
  }

  /// Load on first use, then throw if liblammps is not available.
  void require_loaded() override;

  LammpsLoader(const LammpsLoader &) = delete;
  LammpsLoader &operator=(const LammpsLoader &) = delete;

private:
  LammpsLoader() = default;
  ~LammpsLoader() override;

  void ensure_loaded();

  bool m_loaded{false};
  bool m_tried{false};
  dynlib::Handle m_handle{};
  std::string m_last_error{};

  /// Try to load a symbol; returns nullptr on failure.
  template <typename Fn> Fn load_sym(const char *name) {
    return reinterpret_cast<Fn>(dynlib::sym(m_handle, name));
  }
};

} // namespace eonc
