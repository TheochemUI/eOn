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
#pragma once

#include <memory>

namespace eonc {

class PotRegistry;
class IRAResource;
class ARTnResource;
class PluginLoader;
class MetatomicLoader;

/// Move-only composition root for process resources (dlopen loaders and
/// the potential registry). ClientEON constructs one and Job borrows it;
/// one-shot makeJob owns a unique_ptr<Runtime>. Python Session is this type.
class Runtime {
public:
  Runtime();
  Runtime(Runtime &&) noexcept;
  Runtime &operator=(Runtime &&) noexcept;
  ~Runtime();

  Runtime(const Runtime &) = delete;
  Runtime &operator=(const Runtime &) = delete;

  [[nodiscard]] PotRegistry &pots() noexcept;
  [[nodiscard]] const PotRegistry &pots() const noexcept;
  [[nodiscard]] IRAResource &ira() noexcept;
  [[nodiscard]] const IRAResource &ira() const noexcept;
  [[nodiscard]] ARTnResource &artn() noexcept;
  [[nodiscard]] const ARTnResource &artn() const noexcept;
  [[nodiscard]] PluginLoader &plugins() noexcept;
  [[nodiscard]] const PluginLoader &plugins() const noexcept;
  [[nodiscard]] MetatomicLoader &metatomic() noexcept;
  [[nodiscard]] const MetatomicLoader &metatomic() const noexcept;

private:
  std::unique_ptr<PotRegistry> pots_;
  std::unique_ptr<IRAResource> ira_;
  std::unique_ptr<ARTnResource> artn_;
  std::unique_ptr<PluginLoader> plugins_;
  std::unique_ptr<MetatomicLoader> metatomic_;
};

} // namespace eonc
