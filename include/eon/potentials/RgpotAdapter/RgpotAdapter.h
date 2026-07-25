/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
**
** Adapter over rgpot potential kernels. Migrated pair potentials live in
** librgpot (github.com/OmniPotentRPC/rgpot); eOn's config surface,
** PotRegistry, force counters, and caching stay on this side.
*/
#pragma once

#include <memory>

#include "eon/Potential.h"

#include "rgpot/ForceStructs.hpp"
#include "rgpot/pot_caps.hpp"

/// Wraps one rgpot kernel as an eOn Potential.
///
/// force() passes eOn's flat arrays straight into the kernel's
/// forceImpl — the layouts match field for field, including the
/// row-major 3x3 box — bypassing rgpot's AtomMatrix wrapper and its
/// optional result cache. Thread-sharing policy derives from the
/// kernel's capability descriptor instead of per-type lists.
template <class RPot> class RgpotAdapter final : public Potential {
public:
  /// Kernels construct in place from their config: several hold mutexes
  /// or other immovable state, so the adapter never copies or moves them.
  template <class Cfg>
  RgpotAdapter(PotType ptype, const Parameters &params, const Cfg &cfg)
      : Potential(ptype, params),
        pot_(cfg) {}

  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box) override {
    rgpot::ForceInput in{static_cast<size_t>(N), R, atomicNrs, box};
    rgpot::ForceOut out{F, 0.0, 0.0};
    pot_.forceImpl(in, &out);
    *U = out.energy;
    if (variance != nullptr) {
      *variance = out.variance;
    }
  }

  [[nodiscard]] bool isThreadSafe() const noexcept override {
    return pot_.caps().reentrancy == rgpot::Reentrancy::SharedInstance;
  }

  [[nodiscard]] bool needsPerImageInstance() const noexcept override {
    const auto caps = pot_.caps();
    return caps.perImageInstances ||
           caps.reentrancy == rgpot::Reentrancy::PerInstance;
  }

  [[nodiscard]] const RPot &kernel() const noexcept { return pot_; }

private:
  RPot pot_;
};

/// Factory arm helper: construct the kernel from its config and wrap it.
template <class RPot, class Cfg>
std::shared_ptr<Potential> makeRgpot(PotType ptype, const Parameters &params,
                                     const Cfg &cfg) {
  return std::make_shared<RgpotAdapter<RPot>>(ptype, params, cfg);
}
