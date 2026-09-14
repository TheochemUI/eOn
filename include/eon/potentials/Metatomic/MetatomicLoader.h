/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
**
** Runtime loader for libmetatomic_pot.so (C ABI).
*/
#pragma once

#include "eon/DynLib.h"
#include "metatomic_c_abi.h"

#include <stdexcept>
#include <string>

namespace eonc {

class MetatomicLoader {
public:
  using create_fn = EonMtaPot *(*)(const EonMtaConfig *, char *, size_t);
  using destroy_fn = void (*)(EonMtaPot *);
  using force_fn = int (*)(EonMtaPot *, long, const double *, const int *,
                           double *, double *, double *, const double *);
  using abi_version_fn = int (*)(void);

  static MetatomicLoader &instance();

  /// Attempt (re)load if not yet successful. Safe after potentials_path inject.
  bool try_load();

  [[nodiscard]] bool is_loaded() const noexcept { return m_loaded; }
  void require_loaded();

  create_fn create{nullptr};
  destroy_fn destroy{nullptr};
  force_fn force{nullptr};
  abi_version_fn abi_version{nullptr};

  MetatomicLoader(const MetatomicLoader &) = delete;
  MetatomicLoader &operator=(const MetatomicLoader &) = delete;

private:
  MetatomicLoader() = default;
  ~MetatomicLoader();

  bool m_loaded{false};
  dynlib::Handle m_handle{};
};

} // namespace eonc
