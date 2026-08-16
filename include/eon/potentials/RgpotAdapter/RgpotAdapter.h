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
#include <vector>

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

  /// Kernels with no configuration surface default-construct in place.
  RgpotAdapter(PotType ptype, const Parameters &params)
      : Potential(ptype, params),
        pot_() {}

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

  /// Reports whether the kernel serves a batch natively. rgpot answers
  /// every batch call, falling back to a per-system loop, so this asks
  /// the narrower question the eOn callers care about: is one batched
  /// call cheaper than N single ones? Only a native kernel makes it so.
  [[nodiscard]] bool supportsBatchEvaluation() const noexcept override {
    return pot_.caps().batched;
  }

  /// Hands eOn's batch straight to the kernel. eOn batches images of one
  /// system, so every system shares an atom count; rgpot allows them to
  /// differ, and filling nAtoms per system costs nothing and keeps the
  /// two contracts independent.
  ///
  /// Calls forceBatchImpl rather than forceBatch for the same reason
  /// force() calls forceImpl: eOn owns caching, and routing through
  /// rgpot's cache would key the same geometry twice.
  void forceBatch(long nSystems, long nAtoms, const double *const *positions,
                  const int *const *atomicNrs, double *const *forces,
                  double *energies, double *variances,
                  const double *const *boxes) override {
    const auto n = static_cast<size_t>(nSystems);
    std::vector<rgpot::ForceInput> in;
    std::vector<rgpot::ForceOut> out;
    in.reserve(n);
    out.reserve(n);
    for (size_t i = 0; i < n; ++i) {
      in.push_back(rgpot::ForceInput{static_cast<size_t>(nAtoms), positions[i],
                                     atomicNrs[i], boxes[i]});
      out.push_back(rgpot::ForceOut{forces[i], 0.0, 0.0});
    }

    pot_.forceBatchImpl(
        rgpot::ForceBatch{.nSystems = n, .in = in.data(), .out = out.data()});

    for (size_t i = 0; i < n; ++i) {
      energies[i] = out[i].energy;
      if (variances != nullptr) {
        variances[i] = out[i].variance;
      }
      forceCallCounter++;
      PotRegistry::get().on_force_call(ptype);
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

/// Factory arm helper for kernels whose parameters are fixed tabulated
/// data with no eOn-side configuration surface.
template <class RPot>
std::shared_ptr<Potential> makeRgpotDefault(PotType ptype,
                                            const Parameters &params) {
  return std::make_shared<RgpotAdapter<RPot>>(ptype, params);
}

/// Factory arm helper: construct the kernel from its config and wrap it.
template <class RPot, class Cfg>
std::shared_ptr<Potential> makeRgpot(PotType ptype, const Parameters &params,
                                     const Cfg &cfg) {
  return std::make_shared<RgpotAdapter<RPot>>(ptype, params, cfg);
}
