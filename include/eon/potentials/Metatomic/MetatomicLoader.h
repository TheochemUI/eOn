/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
**
** Runtime loader for libmetatomic_pot.so (C ABI).
*/
#pragma once

#include "eon/DynLib.h"
#include "eon/potentials/Metatomic/metatomic_c_abi.h"

#include <stdexcept>
#include <string>

namespace eonc {

// ---------------------------------------------------------------------------
// IMetatomicLoader: injectable ABI. Production uses instance().
// ---------------------------------------------------------------------------
class IMetatomicLoader {
public:
  using create_fn = EonMtaPot *(*)(const EonMtaConfig *, char *, size_t);
  using destroy_fn = void (*)(EonMtaPot *);
  using force_fn = int (*)(EonMtaPot *, long, const double *, const int *,
                           double *, double *, double *, const double *);
  using abi_version_fn = int (*)(void);

  create_fn create{nullptr};
  destroy_fn destroy{nullptr};
  force_fn force{nullptr};
  abi_version_fn abi_version{nullptr};

  virtual ~IMetatomicLoader() = default;
  virtual void require_loaded() = 0;
  [[nodiscard]] virtual bool is_loaded() const noexcept = 0;

  IMetatomicLoader(const IMetatomicLoader &) = delete;
  IMetatomicLoader &operator=(const IMetatomicLoader &) = delete;

protected:
  IMetatomicLoader() = default;
};

// ---------------------------------------------------------------------------
// MetatomicLoader: process-default loader (Meyer's singleton). Does not
// dlopen until require_loaded() / try_load().
// ---------------------------------------------------------------------------
class MetatomicLoader : public IMetatomicLoader {
public:
  /// Thread-safe singleton accessor (Meyer's pattern). Does not dlopen.
  static MetatomicLoader &instance();

  /// Attempt (re)load if not yet successful. Safe after potentials_path inject.
  bool try_load();

  [[nodiscard]] bool is_loaded() const noexcept override { return m_loaded; }

  /// Load on first use, then throw if libmetatomic_pot is not available.
  void require_loaded() override;

  ~MetatomicLoader() override;
  MetatomicLoader(const MetatomicLoader &) = delete;
  MetatomicLoader &operator=(const MetatomicLoader &) = delete;

private:
  MetatomicLoader() = default;
  friend class Runtime;

  bool m_loaded{false};
  dynlib::Handle m_handle{};
};

} // namespace eonc
