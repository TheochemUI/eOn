/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
**
** Thin host-side Potential: all torch/metatomic work happens in
** libmetatomic_pot.so via C ABI (dlopen).
*/
#pragma once

#include "eon/Potential.h"
#include "eon/potentials/Metatomic/metatomic_c_abi.h"

namespace eonc {

class IMetatomicLoader;

class MetatomicDynPot : public eonc::Potential {
public:
  /// Production: process-default MetatomicLoader.
  explicit MetatomicDynPot(const eonc::Parameters &params);
  /// Test seam: injected loader, no process-default load.
  MetatomicDynPot(const eonc::Parameters &params, IMetatomicLoader &loader);
  ~MetatomicDynPot() override;

  void force(long nAtoms, const double *positions, const int *atomicNrs,
             double *forces, double *energy, double *variance,
             const double *box) override;

  [[nodiscard]] bool isThreadSafe() const noexcept override { return false; }
  [[nodiscard]] bool needsPerImageInstance() const noexcept override {
    return false;
  }

private:
  IMetatomicLoader &loader_;
  EonMtaPot *m_handle{nullptr};
};

} // namespace eonc
