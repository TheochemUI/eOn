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
// IMetatomicLoader: injectable ABI. Production uses
// MetatomicLoader::instance().
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
  /// Attempt (re)load if not yet successful. Safe after potentials_path inject.
  virtual bool try_load() = 0;
  [[nodiscard]] virtual bool is_loaded() const noexcept = 0;
  virtual void require_loaded() = 0;

  IMetatomicLoader(const IMetatomicLoader &) = delete;
  IMetatomicLoader &operator=(const IMetatomicLoader &) = delete;

protected:
  IMetatomicLoader() = default;
};

// ---------------------------------------------------------------------------
// MetatomicLoader: process-default loader (Meyer's singleton). Does not
// dlopen until try_load() / require_loaded().
// ---------------------------------------------------------------------------
class MetatomicLoader : public IMetatomicLoader {
public:
  static MetatomicLoader &instance();

  bool try_load() override;
  [[nodiscard]] bool is_loaded() const noexcept override { return m_loaded; }
  void require_loaded() override;

  MetatomicLoader(const MetatomicLoader &) = delete;
  MetatomicLoader &operator=(const MetatomicLoader &) = delete;

private:
  MetatomicLoader() = default;
  ~MetatomicLoader() override;

  bool m_loaded{false};
  dynlib::Handle m_handle{};
};

} // namespace eonc
