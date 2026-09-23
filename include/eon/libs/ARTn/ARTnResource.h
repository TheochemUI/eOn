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

#include "eon/DynLib.h"
#include <mutex>
#include <stdexcept>
#include <string>

namespace eonc {

// ---------------------------------------------------------------------------
// IARTnResource: injectable ABI. Production uses ARTnResource::instance().
// ---------------------------------------------------------------------------
class IARTnResource {
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
  /// pARTn bind(C) get_param / get_runparam: type(c_ptr), intent(out).
  /// That is void** in C. The artn.h comment says void*. A void* out
  /// parameter writes the allocated pointer at the wrong address.
  using get_param_fn = int (*)(const char *name, void **cval);
  using get_runparam_fn = int (*)(const char *name, void **cval);
  using get_data_fn = int (*)(const char *name, void **cval);
  using print_caller_fn = void (*)();
  using artn_step_fn = void (*)(const int nat, const double etot,
                                double *const force, int const *ityp,
                                double *const pos, const double *box,
                                const int *if_pos, double *displ_vec,
                                bool *lconv);
  /// Retrieve pARTn error code + c_malloc'd message pointer.
  /// Defined in m_artn_error.f90 as bind(C, name="get_error"):
  ///   int get_error(void** cmsg);
  /// Non-zero return indicates an error; *cmsg then points at a
  /// C-string the caller must std::free. Zero leaves *cmsg null.
  using get_error_fn = int (*)(void **cmsg);

  /// Serializes access to the Fortran backend's shared global state.
  std::mutex library_mutex;

  virtual ~IARTnResource() = default;
  virtual void require_loaded() = 0;
  [[nodiscard]] virtual bool is_loaded() const noexcept = 0;
  [[nodiscard]] virtual artn_create_fn get_create_fn() const = 0;
  [[nodiscard]] virtual setup_artn_fn get_setup_fn() const = 0;
  [[nodiscard]] virtual artn_fn get_artn_fn() const = 0;
  [[nodiscard]] virtual artn_destroy_fn get_destroy_fn() const = 0;
  [[nodiscard]] virtual set_param_fn get_set_param_fn() const = 0;
  [[nodiscard]] virtual get_param_fn get_get_param_fn() const = 0;
  [[nodiscard]] virtual get_runparam_fn get_get_runparam_fn() const = 0;
  [[nodiscard]] virtual get_data_fn get_get_data_fn() const = 0;
  [[nodiscard]] virtual print_caller_fn get_print_caller_fn() const = 0;
  [[nodiscard]] virtual artn_step_fn get_artn_step_fn() const = 0;
  /// May be null on older pARTn builds that predate the C get_error
  /// wrapper; callers must null-check before dispatching.
  [[nodiscard]] virtual get_error_fn get_get_error_fn() const = 0;

  IARTnResource(const IARTnResource &) = delete;
  IARTnResource &operator=(const IARTnResource &) = delete;

protected:
  IARTnResource() = default;
};

// ---------------------------------------------------------------------------
// ARTnResource: process-default loader (Meyer's singleton).
// ---------------------------------------------------------------------------
class ARTnResource : public IARTnResource {
public:
  /// Thread-safe singleton accessor (Meyer's pattern).
  static ARTnResource &instance();

  [[nodiscard]] bool is_loaded() const noexcept override { return m_loaded; }

  /// Throws std::runtime_error if libartn is not available.
  void require_loaded() override;

  [[nodiscard]] artn_create_fn get_create_fn() const override {
    return artn_create_;
  }
  [[nodiscard]] setup_artn_fn get_setup_fn() const override {
    return setup_artn_;
  }
  [[nodiscard]] artn_fn get_artn_fn() const override { return artn_; }
  [[nodiscard]] artn_destroy_fn get_destroy_fn() const override {
    return artn_destroy_;
  }
  [[nodiscard]] set_param_fn get_set_param_fn() const override {
    return set_param_;
  }
  [[nodiscard]] get_param_fn get_get_param_fn() const override {
    return get_param_;
  }
  [[nodiscard]] get_runparam_fn get_get_runparam_fn() const override {
    return get_runparam_;
  }
  [[nodiscard]] get_data_fn get_get_data_fn() const override {
    return get_data_;
  }
  [[nodiscard]] print_caller_fn get_print_caller_fn() const override {
    return print_caller_;
  }
  [[nodiscard]] artn_step_fn get_artn_step_fn() const override {
    return artn_step_;
  }
  [[nodiscard]] get_error_fn get_get_error_fn() const override {
    return get_error_;
  }

  ~ARTnResource() override;
  ARTnResource(const ARTnResource &) = delete;
  ARTnResource &operator=(const ARTnResource &) = delete;

private:
  ARTnResource();
  friend class Runtime;

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
  get_error_fn get_error_{nullptr};

  /// Try to load a symbol; returns nullptr on failure.
  template <typename Fn> Fn load_sym(const char *name) const {
    return reinterpret_cast<Fn>(dynlib::sym(m_handle, name));
  }
};

/// Global access to thread-safe ARTn resource
inline ARTnResource &get_artn_resource() { return ARTnResource::instance(); }

} // namespace eonc
