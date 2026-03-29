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
#include "EonLogger.h"

#include "Eigen.h"
#include "Parameters.h"
#include "PotRegistry.h"
#include <atomic>
#include <memory>

namespace eonc {

class Potential {
protected:
  PotType ptype;

private:
  uint64_t m_registry_id;
  PotRegistry::TimePoint m_created_at;

public:
  std::atomic<size_t> forceCallCounter;

  // Main Constructor (no Parameters dependency)
  explicit Potential(PotType a_ptype)
      : ptype{a_ptype},
        m_registry_id{PotRegistry::get().on_created(a_ptype)},
        m_created_at{PotRegistry::Clock::now()},
        forceCallCounter{0} {}

  // Convenience constructor from Parameters (for backward compat)
  Potential(PotType a_ptype, const Parameters &)
      : Potential(a_ptype) {}

  Potential(const Parameters &a_params)
      : Potential(a_params.potential_options.potential) {}

  virtual ~Potential() {
    PotRegistry::get().on_destroyed(m_registry_id, ptype, forceCallCounter,
                                    m_created_at);
  }

  // Does not take into account the fixed / free atoms
  // Variance here is null when not needed and that's OK
  void virtual force(long nAtoms, const double *positions, const int *atomicNrs,
                     double *forces, double *energy, double *variance,
                     const double *box) = 0;

  std::tuple<double, AtomMatrix>
  get_ef(const AtomMatrix &pos, const VectorXi &atmnrs, const Matrix3d &box);

  [[nodiscard]] PotType getType() const { return this->ptype; }

  /// Whether this is a surrogate (GP) potential. Override in
  /// SurrogatePotential.
  [[nodiscard]] virtual bool isSurrogate() const noexcept { return false; }

  /// Whether this potential's force() can be called from multiple threads
  /// on the SAME instance. Python-based potentials return false.
  /// Potentials with internal mutex (MetatomicPotential) return true but
  /// serialize internally -- use needsPerImageInstance() to check if
  /// separate instances would enable true parallelism.
  [[nodiscard]] virtual bool isThreadSafe() const noexcept { return true; }

  /// Whether NEB should create separate Potential instances per image
  /// for true parallel force evaluation. When true, NEB calls
  /// makePotential() once per image instead of sharing one instance.
  /// Override in potentials that use internal mutexes (e.g.
  /// MetatomicPotential).
  [[nodiscard]] virtual bool needsPerImageInstance() const noexcept {
    return false;
  }

  /// Whether this potential supports batched evaluation of N systems in a
  /// single call. When true, callers (NEB, Dimer) should use forceBatch()
  /// instead of N individual force() calls for better GPU utilization.
  [[nodiscard]] virtual bool supportsBatchEvaluation() const noexcept {
    return false;
  }

  /// Evaluate forces for N systems in a single call. Default: loops over
  /// force(). Override in potentials that support native batching (e.g.
  /// MetatomicPotential uses a single model.forward() for all N systems).
  virtual void forceBatch(long nSystems, long nAtoms,
                          const double *const *positions,
                          const int *const *atomicNrs, double *const *forces,
                          double *energies, double *variances,
                          const double *const *boxes) {
    for (long i = 0; i < nSystems; i++) {
      double var = 0;
      force(nAtoms, positions[i], atomicNrs[i], forces[i], &energies[i], &var,
            boxes[i]);
      if (variances)
        variances[i] = var;
      forceCallCounter++;
      PotRegistry::get().on_force_call(ptype);
    }
  }
};

namespace helpers {
std::shared_ptr<Potential> makePotential(const Parameters &params);
std::shared_ptr<Potential> makePotential(PotType ptype,
                                         const Parameters &params);
} // namespace helpers

} // namespace eonc

using eonc::Potential;
