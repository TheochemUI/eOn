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

/// RAII resource manager for the ARTn C library with global synchronization.
///
/// Uses dlopen (POSIX) or LoadLibrary (Windows) to load libartn at runtime.
/// Implements global mutex to ensure thread-safety when multiple threads
/// access the Fortran backend which has shared global state.
///
/// WARNING: Due to the non-thread-safe Fortran backend, all ARTn operations
/// are serialized through a global mutex. Multiple threads calling ARTn
/// simultaneously will execute sequentially, potentially creating a significant
/// performance bottleneck in multi-threaded contexts. Consider process-level
/// parallelism if true concurrent ARTn operations are required.

#include "../../DynLib.h"
#include <mutex>
#include <stdexcept>
#include <string>

namespace eonc {

class ARTnResource {
public:
  // ARTn C API function pointer types (from artn.h)
  using artn_create_fn = int (*)();
  using setup_artn_fn = void (*)(const int nat, bool *cerr);
  using artn_fn = void (*)(const int nat, const double *etot, const double *f,
                           int const *ityp, double *const tau, const int *order,
                           const double *lat, const int *if_pos, int *disp_code,
                           double *disp_vec, bool *lconv);
  using artn_destroy_fn = void (*)();
  using set_param_fn = int (*)(const char *const name, const int crank,
                               const int *csize, const void *cval);
  using get_param_fn = int (*)(const char *name, void *cval);
  using get_runparam_fn = int (*)(const char *name, void *cval);
  using get_data_fn = int (*)(const char *name, void **cval);
  using print_caller_fn = void (*)();
  using artn_step_fn = void (*)(const int nat, const double etot,
                                double *const force, int const *ityp,
                                double *const pos, const double *box,
                                const int *if_pos, double *displ_vec,
                                bool *lconv);

  /// Singleton accessor (Meyer's pattern).
  static ARTnResource &instance();

  /// Global mutex to ensure only one thread accesses ARTn library at a time
  std::mutex library_mutex;

  /// True if libartn was successfully loaded.
  [[nodiscard]] bool is_loaded() const noexcept { return m_loaded; }

  /// Throws std::runtime_error if libartn is not available.
  void require_loaded() const;

  // Function pointer accessors
  artn_create_fn get_create_fn() const { return artn_create_; }
  setup_artn_fn get_setup_fn() const { return setup_artn_; }
  artn_fn get_artn_fn() const { return artn_; }
  artn_destroy_fn get_destroy_fn() const { return artn_destroy_; }
  set_param_fn get_set_param_fn() const { return set_param_; }
  get_param_fn get_get_param_fn() const { return get_param_; }
  get_runparam_fn get_get_runparam_fn() const { return get_runparam_; }
  get_data_fn get_get_data_fn() const { return get_data_; }
  print_caller_fn get_print_caller_fn() const { return print_caller_; }
  artn_step_fn get_artn_step_fn() const { return artn_step_; }

  /// Non-copyable, non-movable
  ARTnResource(const ARTnResource &) = delete;
  ARTnResource &operator=(const ARTnResource &) = delete;

private:
  ARTnResource();
  ~ARTnResource();

  bool m_loaded{false};
  dynlib::Handle m_handle{};

  // Loaded function pointers (null if library not found)
  artn_create_fn artn_create_{nullptr};
  setup_artn_fn setup_artn_{nullptr};
  artn_fn artn_{nullptr};
  artn_destroy_fn artn_destroy_{nullptr};
  set_param_fn set_param_{nullptr};
  get_param_fn get_param_{nullptr};
  get_runparam_fn get_runparam_{nullptr};
  get_data_fn get_data_{nullptr};
  print_caller_fn print_caller_{nullptr};
  artn_step_fn artn_step_{nullptr};

  /// Try to load a symbol; returns nullptr on failure.
  template <typename Fn> Fn load_sym(const char *name) const {
    return reinterpret_cast<Fn>(dynlib::sym(m_handle, name));
  }
};

/// Global access to thread-safe ARTn resource
inline ARTnResource &get_artn_resource() { return ARTnResource::instance(); }

} // namespace eonc
