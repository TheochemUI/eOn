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

#include "ARTnResource.h"
#include <iostream>

namespace eonc {

ARTnResource &ARTnResource::instance() {
  static ARTnResource resource;
  return resource;
}

ARTnResource::ARTnResource() {
#ifdef _WIN32
  const char *names[] = {"artn.dll", "libartn.dll", nullptr};
#elif defined(__APPLE__)
  const char *names[] = {"libartn.dylib", "libartn.0.dylib", nullptr};
#else
  const char *names[] = {"libartn.so", "libartn.so.0", nullptr};
#endif

  m_handle = dynlib::openFirst(names);
  if (!m_handle) {
    return; // Not found -- require_loaded() will throw if user requests ARTn
  }

  // Load required symbols. eOn's lifecycle uses artn_destroy for teardown;
  // clean_artn (the old name kept in pARTn) is intentionally not loaded here.
  artn_create_ = load_sym<artn_create_fn>("artn_create");
  setup_artn_ = load_sym<setup_artn_fn>("setup_artn");
  artn_ = load_sym<artn_fn>("artn");
  artn_destroy_ = load_sym<artn_destroy_fn>("artn_destroy");
  set_param_ = load_sym<set_param_fn>("set_param");
  get_param_ = load_sym<get_param_fn>("get_param");
  get_runparam_ = load_sym<get_runparam_fn>("get_runparam");
  get_data_ = load_sym<get_data_fn>("get_data");
  print_caller_ = load_sym<print_caller_fn>("print_caller");
  artn_step_ = load_sym<artn_step_fn>("artn_step");
  // Optional: older artn-plugin builds (pre-get_error C wrapper) may not
  // export this. Callers null-check via get_get_error_fn() before use.
  get_error_ = load_sym<get_error_fn>("get_error");

  // Check that all required symbols are loaded
  if (!artn_create_ || !setup_artn_ || !artn_ || !artn_destroy_ ||
      !set_param_ || !get_param_ || !get_runparam_ || !get_data_ ||
      !print_caller_ || !artn_step_) {
    // std::cerr instead of quill: this runs from the static-initialiser of
    // instance() before eOn's logger backend is guaranteed to be configured.
    std::cerr << "[ARTN] Library loaded but missing required symbols\n";
    dynlib::close(m_handle);
    m_handle = {};
    return;
  }

  m_loaded = true;
}

ARTnResource::~ARTnResource() {
  // Intentionally do NOT dlclose at static destruction. libartn pulls in
  // gfortran runtime + LAPACK; its destructor chain runs Fortran
  // FINAL/STOP statements that touch stdio that is already being torn
  // down by the time this Meyer-singleton dtor runs. The handle is
  // process-lifetime; OS unmaps the lib at exit.
}

void ARTnResource::require_loaded() const {
  if (!m_loaded) {
    throw std::runtime_error(
        "ARTn saddle search requested but libartn not found.\n"
        "Install the ARTn library and ensure it's in your library search path "
        "(LD_LIBRARY_PATH / DYLD_LIBRARY_PATH / PATH).");
  }
}

} // namespace eonc
