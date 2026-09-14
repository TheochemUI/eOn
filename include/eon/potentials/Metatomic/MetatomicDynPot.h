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
#include "metatomic_c_abi.h"

namespace eonc {

class MetatomicDynPot : public Potential {
public:
  explicit MetatomicDynPot(const Parameters &params);
  ~MetatomicDynPot() override;

  void force(long nAtoms, const double *positions, const int *atomicNrs,
             double *forces, double *energy, double *variance,
             const double *box) override;

  [[nodiscard]] bool isThreadSafe() const noexcept override { return false; }
  [[nodiscard]] bool needsPerImageInstance() const noexcept override {
    return false;
  }

private:
  EonMtaPot *m_handle{nullptr};
};

} // namespace eonc
