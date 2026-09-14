/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
*/
#include "eon/potentials/Metatomic/MetatomicLoader.h"

#include "eon/potentials/PluginLoader.h"

#include <iostream>
#include <vector>

#ifndef _WIN32
#include <dlfcn.h>
#endif

namespace eonc {

MetatomicLoader &MetatomicLoader::instance() {
  static MetatomicLoader loader;
  return loader;
}

bool MetatomicLoader::try_load() {
  if (m_loaded)
    return true;

#ifdef _WIN32
  const char *candidates[] = {"metatomic_pot.dll", "libmetatomic_pot.dll",
                              nullptr};
#elif defined(__APPLE__)
  const char *candidates[] = {"libmetatomic_pot.dylib", "libmetatomic_pot.so",
                              nullptr};
#else
  const char *candidates[] = {"libmetatomic_pot.so", nullptr};
#endif

  // RTLD_GLOBAL without DEEPBIND so MetatomicPotential's Potential/PotRegistry
  // symbols bind to the already-loaded host libeonclib (single registry).
  auto open_mta = [](const char *name) -> dynlib::Handle {
#if defined(__linux__)
    return ::dlopen(name, RTLD_NOW | RTLD_GLOBAL);
#elif defined(_WIN32)
    return dynlib::open(name);
#else
    return ::dlopen(name, RTLD_NOW | RTLD_GLOBAL);
#endif
  };

  for (const char **p = candidates; *p; ++p) {
    m_handle = open_mta(*p);
    if (m_handle)
      break;
  }
  if (!m_handle) {
    for (const auto &dir : PluginLoader::instance().search_paths()) {
      for (const char **p = candidates; *p; ++p) {
        std::string full = dir;
        if (!full.empty() && full.back() != '/' && full.back() != '\\')
          full += '/';
        full += *p;
        m_handle = open_mta(full.c_str());
        if (m_handle)
          break;
      }
      if (m_handle)
        break;
    }
  }
  if (!m_handle)
    return false;

  create = dynlib::loadSym<create_fn>(m_handle, "eon_mta_pot_create");
  destroy = dynlib::loadSym<destroy_fn>(m_handle, "eon_mta_pot_destroy");
  force = dynlib::loadSym<force_fn>(m_handle, "eon_mta_pot_force");
  abi_version =
      dynlib::loadSym<abi_version_fn>(m_handle, "eon_mta_abi_version");

  if (!create || !destroy || !force || !abi_version) {
    std::cerr << "[Metatomic] libmetatomic_pot loaded but missing C ABI "
                 "symbols\n";
    dynlib::close(m_handle);
    m_handle = {};
    create = nullptr;
    destroy = nullptr;
    force = nullptr;
    abi_version = nullptr;
    return false;
  }
  const int ver = abi_version();
  if (ver != EON_MTA_ABI_VERSION) {
    std::cerr << "[Metatomic] ABI version mismatch: library=" << ver
              << " host=" << EON_MTA_ABI_VERSION << "\n";
    dynlib::close(m_handle);
    m_handle = {};
    create = nullptr;
    destroy = nullptr;
    force = nullptr;
    abi_version = nullptr;
    return false;
  }
  m_loaded = true;
  return true;
}

MetatomicLoader::~MetatomicLoader() { dynlib::close(m_handle); }

void MetatomicLoader::require_loaded() {
  if (!try_load()) {
    throw std::runtime_error(
        "Metatomic potential requested but libmetatomic_pot.so not found or "
        "ABI mismatch.\n"
        "Build with -Dwith_metatomic=true and put the plugin on "
        "LD_LIBRARY_PATH / EON_POTENTIALS_PATH / [Potential] potentials_path.");
  }
}

} // namespace eonc
