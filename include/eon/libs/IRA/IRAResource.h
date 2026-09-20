/*
 * This file is part of eOn.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 *
 * Copyright (c) 2010--present, eOn Development Team
 * All rights reserved.
 *
 * Repo:
 * https://github.com/TheochemUI/eOn
 */

#pragma once

/// RAII resource manager for the IRA C library with global synchronization.
///
/// Uses dlopen (POSIX) or LoadLibrary (Windows) to load libira at runtime.
/// Implements global mutex to ensure thread-safety when multiple threads
/// access the Fortran backend which has shared global state.
///
/// WARNING: Due to the non-thread-safe Fortran backend, all IRA operations
/// are serialized through a global mutex. Multiple threads calling IRA
/// simultaneously will execute sequentially, potentially creating a significant
/// performance bottleneck in multi-threaded contexts. Consider process-level
/// parallelism if true concurrent IRA operations are required.

#include "eon/DynLib.h"
#include <mutex>
#include <stdexcept>
#include <string>

namespace eonc {

// ---------------------------------------------------------------------------
// IIRAResource: injectable ABI. Production uses IRAResource::instance().
// ---------------------------------------------------------------------------
class IIRAResource {
public:
  // IRA C API function pointer types (from iralib_interf.h)
  using libira_match_fn = void (*)(int nat1, const int *typ1,
                                   const double *coords1, const int *cand1,
                                   int nat2, const int *typ2,
                                   const double *coords2, const int *cand2,
                                   double distThreshold, double **rotMat,
                                   double **trans, int **perm,
                                   double *hausdorffDist, int *ierr);
  using libira_cshda_pbc_fn = void (*)(int nat1, const int *typ1,
                                       const double *coords1, int nat2,
                                       const int *typ2, const double *coords2,
                                       const double *lat, double distThreshold,
                                       int **found, double **dists);
  using libira_compute_all_fn =
      void (*)(int nat, const int *typ, const double *coords, double threshold,
               int prescreenIh, int *n_mat, double **mat, int **perm, char **op,
               int **n, int **p, double **ax, double **angle, double **dH,
               char **pg, int *n_prin_ax, double **prin_ax, int *cerr);
  using libira_get_nmax_fn = int (*)();

  /// Serializes access to the Fortran backend's shared global state.
  std::mutex library_mutex;

  virtual ~IIRAResource() = default;
  virtual void require_loaded() = 0;
  [[nodiscard]] virtual bool is_loaded() const noexcept = 0;
  [[nodiscard]] virtual libira_match_fn get_match_fn() const = 0;
  [[nodiscard]] virtual libira_cshda_pbc_fn get_cshda_pbc_fn() const = 0;
  [[nodiscard]] virtual libira_compute_all_fn get_compute_all_fn() const = 0;
  [[nodiscard]] virtual libira_get_nmax_fn get_get_nmax_fn() const = 0;

  IIRAResource(const IIRAResource &) = delete;
  IIRAResource &operator=(const IIRAResource &) = delete;

protected:
  IIRAResource() = default;
};

// ---------------------------------------------------------------------------
// IRAResource: process-default loader (Meyer's singleton).
// ---------------------------------------------------------------------------
class IRAResource : public IIRAResource {
public:
  /// Thread-safe singleton accessor (Meyer's pattern).
  static IRAResource &instance();

  [[nodiscard]] bool is_loaded() const noexcept override { return m_loaded; }

  /// Throws std::runtime_error if libira is not available.
  void require_loaded() override;

  [[nodiscard]] libira_match_fn get_match_fn() const override {
    return libira_match_;
  }
  [[nodiscard]] libira_cshda_pbc_fn get_cshda_pbc_fn() const override {
    return libira_cshda_pbc_;
  }
  [[nodiscard]] libira_compute_all_fn get_compute_all_fn() const override {
    return libira_compute_all_;
  }
  [[nodiscard]] libira_get_nmax_fn get_get_nmax_fn() const override {
    return libira_get_nmax_;
  }

  ~IRAResource() override;
  IRAResource(const IRAResource &) = delete;
  IRAResource &operator=(const IRAResource &) = delete;

private:
  IRAResource();
  friend class Runtime;

  bool m_loaded{false};
  dynlib::Handle m_handle{};

  // Loaded function pointers (null if library not found)
  libira_match_fn libira_match_{nullptr};
  libira_cshda_pbc_fn libira_cshda_pbc_{nullptr};
  libira_compute_all_fn libira_compute_all_{nullptr};
  libira_get_nmax_fn libira_get_nmax_{nullptr};

  /// Try to load a symbol; returns nullptr on failure.
  template <typename Fn> Fn load_sym(const char *name) const {
    return reinterpret_cast<Fn>(dynlib::sym(m_handle, name));
  }
};

/// Global access to thread-safe IRA resource
inline IRAResource &get_ira_resource() { return IRAResource::instance(); }

} // namespace eonc
