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

#include "Optimizer.h"
#include "client/matter/Matter.h"

namespace eonc {
class MinimizationJob {
  struct Params final {
    // We don't know what comes here. It is a bad idea to have this enumerated
    // at compile time, since the idea is that the optimizers can be extended at
    // runtime, e.g. via Python, so no std::variant tricks with monostate either
    // It is always possible to instead give up and use TOML as the defacto way
    // to handle any kind of optimizer, but that is unpleasant as well as would
    // require a stronger dependence on tomlpp than we have now.
    //
    // Would prevent using any other kind of input format.
    //
    // It is still probably OK for now, to try to use a variant for the
    // optimizers, but there will be trouble later probably.
    //
    // On the other hand, using variant for the optimizers will enable handling
    // optimizers as a first class object again (value semantics)
    //
    // Recall that the issue here is that for "final" parameter packs like
    // Potential, things are OK as is, and even for nested parameter packs like
    // Optimizer.CG
    //
    // The problem and question is what happens for Jobs like MinimizationJob
    // which are required to know both their own ([Main]) parameters, and "pass
    // through" certain sections of the configuration without either:
    //
    // - Hard dependency on TOML
    // - Too many layers of indirection (i.e. MinConfig struct proxying)
    // - std::variant on the Parameters of the optimizer, as an example
    //
    // The last case is bad because there's no guarantee that there will only be
    // two, what happens when we want to additionally pass say, a [Debug]
    // section?
    //
    // These are all stemming from the fact that unlike [Potential] which never
    // constructed beyond what we know upfront, the Jobs are free to construct
    // sub-objects, which need their parameters, preferrably without much
    // pointer hacking or storing another Parameters god object...
  };

public:
  MinimizationJob(const OptimBase::Params &_p)
      : m_p{_p} {
    m_log = spdlog::get("combi");
  }
  ~MinimizationJob(void) = default;
  bool runImpl(Matter &);

private:
  OptimBase::Params m_p;
  std::shared_ptr<spdlog::logger> m_log;
};

} // namespace eonc
