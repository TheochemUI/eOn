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
#include "eon/Runtime.h"

#include "eon/PotRegistry.h"
#include "eon/libs/ARTn/ARTnResource.h"
#include "eon/libs/IRA/IRAResource.h"
#include "eon/potentials/Metatomic/MetatomicLoader.h"
#include "eon/potentials/PluginLoader.h"

namespace eonc {

Runtime::Runtime()
    : pots_(std::make_unique<PotRegistry>()),
      ira_(std::unique_ptr<IRAResource>(new IRAResource())),
      artn_(std::unique_ptr<ARTnResource>(new ARTnResource())),
      plugins_(std::unique_ptr<PluginLoader>(new PluginLoader())),
      metatomic_(std::unique_ptr<MetatomicLoader>(new MetatomicLoader())) {}

Runtime::Runtime(Runtime &&) noexcept = default;
Runtime &Runtime::operator=(Runtime &&) noexcept = default;
Runtime::~Runtime() = default;

PotRegistry &Runtime::pots() noexcept { return *pots_; }
const PotRegistry &Runtime::pots() const noexcept { return *pots_; }
IRAResource &Runtime::ira() noexcept { return *ira_; }
const IRAResource &Runtime::ira() const noexcept { return *ira_; }
ARTnResource &Runtime::artn() noexcept { return *artn_; }
const ARTnResource &Runtime::artn() const noexcept { return *artn_; }
PluginLoader &Runtime::plugins() noexcept { return *plugins_; }
const PluginLoader &Runtime::plugins() const noexcept { return *plugins_; }
MetatomicLoader &Runtime::metatomic() noexcept { return *metatomic_; }
const MetatomicLoader &Runtime::metatomic() const noexcept {
  return *metatomic_;
}

} // namespace eonc
