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

#include "IRAResource.h"
#include <iostream>

namespace eonc {

IRAResource &IRAResource::instance() {
  static IRAResource resource;
  return resource;
}

IRAResource::IRAResource() {
#ifdef _WIN32
  const char *names[] = {"ira.dll", "libira.dll", nullptr};
#elif defined(__APPLE__)
  const char *names[] = {"libira.dylib", "libira.0.dylib", nullptr};
#else
  const char *names[] = {"libira.so", "libira.so.0", nullptr};
#endif

  m_handle = dynlib::openFirst(names);
  if (!m_handle) {
    return; // Not found -- require_loaded() will throw if user requests IRA
  }

  // Load required symbols
  libira_match_ = load_sym<libira_match_fn>("libira_match");
  libira_cshda_pbc_ = load_sym<libira_cshda_pbc_fn>("libira_cshda_pbc");
  libira_compute_all_ = load_sym<libira_compute_all_fn>("libira_compute_all");
  libira_get_nmax_ = load_sym<libira_get_nmax_fn>("libira_get_nmax");

  // Check that all required symbols are loaded
  if (!libira_match_ || !libira_cshda_pbc_ || !libira_compute_all_ ||
      !libira_get_nmax_) {
    // std::cerr instead of quill: this runs from the static-initialiser of
    // instance() before eOn's logger backend is guaranteed to be configured.
    std::cerr << "[IRA] Library loaded but missing required symbols\n";
    dynlib::close(m_handle);
    m_handle = {};
    return;
  }

  m_loaded = true;
}

IRAResource::~IRAResource() {
  // Intentionally do NOT dlclose at static destruction. libira is
  // Fortran with its own runtime fini chain; mirroring the XtbLoader /
  // LammpsLoader / ARTnResource policy of letting the OS reclaim the
  // mapping at process exit.
}

void IRAResource::require_loaded() const {
  if (!m_loaded) {
    throw std::runtime_error(
        "IRA structure comparison requested but libira not found.\n"
        "Install the IRA library and ensure it's in your library search path "
        "(LD_LIBRARY_PATH / DYLD_LIBRARY_PATH / PATH).");
  }
}

} // namespace eonc
