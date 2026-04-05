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

  // Load required symbols
  artn_create_ = load_sym<artn_create_fn>("artn_create");
  setup_artn_ = load_sym<setup_artn_fn>("setup_artn");
  artn_ = load_sym<artn_fn>("artn");
  clean_artn_ = load_sym<clean_artn_fn>("clean_artn");
  set_param_ = load_sym<set_param_fn>("set_param");
  get_param_ = load_sym<get_param_fn>("get_param");
  get_runparam_ = load_sym<get_runparam_fn>("get_runparam");
  get_data_ = load_sym<get_data_fn>("get_data");
  print_caller_ = load_sym<print_caller_fn>("print_caller");
  artn_step_ = load_sym<artn_step_fn>("artn_step");

  // Check that all required symbols are loaded
  if (!artn_create_ || !setup_artn_ || !artn_ || !clean_artn_ || !set_param_ ||
      !get_param_ || !get_runparam_ || !get_data_ || !print_caller_ ||
      !artn_step_) {
    std::cerr << "[ARTN] Library loaded but missing required symbols\n";
    dynlib::close(m_handle);
    m_handle = {};
    return;
  }

  m_loaded = true;
}

ARTnResource::~ARTnResource() {
  if (m_handle) {
    dynlib::close(m_handle);
    m_handle = {};
  }
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
