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
#include "eon/potentials/Rgpot/GenericEngineLoader.h"
#include <cstdlib>
#include <stdexcept>
#include <vector>
#ifndef _WIN32
#include <dlfcn.h>
#endif

namespace {
void *open_lib(const char *path) {
#ifndef _WIN32
  return dlopen(path, RTLD_NOW | RTLD_GLOBAL);
#else
  return nullptr;
#endif
}
void close_lib(void *h) {
#ifndef _WIN32
  if (h)
    dlclose(h);
#endif
}
void *load_sym(void *h, const char *n) {
#ifndef _WIN32
  return dlsym(h, n);
#else
  return nullptr;
#endif
}
} // namespace

GenericEngineLoader::GenericEngineLoader(const GenericEngineOptions &opt)
    : m_tag(opt.tag.empty() ? opt.library : opt.tag) {
  std::vector<std::string> paths;
  if (!opt.engine_path.empty())
    paths.push_back(opt.engine_path);
  if (!opt.env_var.empty())
    if (const char *e = std::getenv(opt.env_var.c_str()))
      if (e && *e)
        paths.emplace_back(e);
  paths.push_back(opt.library);
  auto add_dirs = [&](const char *env) {
    if (!env)
      return;
    std::string s(env);
    size_t i = 0;
    while (i < s.size()) {
      const auto j = s.find(':', i);
      const auto dir = s.substr(i, j == std::string::npos ? j : j - i);
      if (!dir.empty())
        paths.push_back(dir + "/" + opt.library);
      if (j == std::string::npos)
        break;
      i = j + 1;
    }
  };
  add_dirs(std::getenv("EON_POTENTIALS_PATH"));

  std::string last_dlerr;
  for (const auto &p : paths) {
    m_lib = open_lib(p.c_str());
    if (m_lib)
      break;
#ifndef _WIN32
    if (const char *e = dlerror())
      last_dlerr = e;
#endif
  }
  if (!m_lib) {
    std::string msg =
        "RGPOT(" + m_tag + "): " + opt.library + " not found (set " +
        (opt.env_var.empty() ? std::string("engine_path") : opt.env_var) +
        " or EON_POTENTIALS_PATH)";
    if (!last_dlerr.empty())
      msg += std::string("; last dlerror: ") + last_dlerr;
    throw std::runtime_error(msg);
  }

  auto abi =
      reinterpret_cast<int (*)()>(load_sym(m_lib, "rgpot_engine_abi_version"));
  m_create =
      reinterpret_cast<create_fn>(load_sym(m_lib, "rgpot_engine_create"));
  m_destroy =
      reinterpret_cast<destroy_fn>(load_sym(m_lib, "rgpot_engine_destroy"));
  m_force = reinterpret_cast<force_fn>(load_sym(m_lib, "rgpot_engine_force"));
  if (!abi || !m_create || !m_destroy || !m_force ||
      abi() != RGPOT_ENGINE_ABI_VERSION) {
    close_lib(m_lib);
    m_lib = nullptr;
    throw std::runtime_error("RGPOT(" + m_tag +
                             "): engine C ABI missing/mismatch");
  }
  m_set_cache =
      reinterpret_cast<set_cache_fn>(load_sym(m_lib, "rgpot_engine_set_cache"));
  m_cache_counts = reinterpret_cast<cache_counts_fn>(
      load_sym(m_lib, "rgpot_engine_cache_counts"));
  char err[1024]{};
  m_pot = m_create(opt.config.data(), opt.config.size(), err, sizeof err);
  if (!m_pot) {
    close_lib(m_lib);
    m_lib = nullptr;
    throw std::runtime_error("RGPOT(" + m_tag +
                             "): create failed: " + std::string(err));
  }
}

GenericEngineLoader::~GenericEngineLoader() {
  // Destroy the pot while the engine is still mapped. Do not dlclose the
  // engine: torch static teardown after dlclose routinely SEGV on exit.
  if (m_pot && m_destroy)
    m_destroy(m_pot);
  m_pot = nullptr;
  m_lib = nullptr;
}

std::string GenericEngineLoader::lastError(const char *operation) const {
  using last_error_fn = size_t (*)(const RgpotEnginePot *, char *, size_t);
  auto last = reinterpret_cast<last_error_fn>(
      load_sym(m_lib, "rgpot_engine_last_error"));
  char detail[1024]{};
  if (last)
    last(m_pot, detail, sizeof detail);
  return "RGPOT(" + m_tag + "): " + operation +
         (detail[0] ? ": " + std::string(detail) : "");
}

std::string GenericEngineLoader::version() const {
  const auto version = reinterpret_cast<const char *(*)()>(
      load_sym(m_lib, "rgpot_engine_version"));
  const char *value = version ? version() : nullptr;
  return value ? value : "";
}

void GenericEngineLoader::attachCache(const std::string &path) {
  if (!m_set_cache || !m_cache_counts)
    throw std::runtime_error("RGPOT(" + m_tag +
                             "): engine does not export the result cache ABI");
  if (m_set_cache(m_pot, path.c_str()) != 0)
    throw std::runtime_error(lastError("attach result cache"));
}

void GenericEngineLoader::detachCache() noexcept {
  if (m_set_cache)
    m_set_cache(m_pot, nullptr);
}

std::array<std::uint64_t, 2> GenericEngineLoader::cacheCounts() const {
  std::array<std::uint64_t, 2> counts{};
  if (!m_cache_counts || m_cache_counts(m_pot, &counts[0], &counts[1]) != 0)
    throw std::runtime_error(lastError("read result cache counters"));
  return counts;
}

void GenericEngineLoader::force(long N, const double *R, const int *atomicNrs,
                                double *F, double *U, double *variance,
                                const double *box) const {
  if (!m_pot || !m_force)
    throw std::runtime_error("RGPOT(" + m_tag + "): engine not available");
  double var = 0.0;
  const int rc = m_force(m_pot, N, R, atomicNrs, F, U,
                         variance ? variance : &var, box, nullptr, nullptr);
  if (rc != 0)
    throw std::runtime_error(lastError("force failed"));
}
