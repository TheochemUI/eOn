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

#include "eon/Matter.h"
#include "eon/Potential.h"
#include "eon/relax/eon_relax_engine.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <vector>

namespace eonc {

struct SurfaceRecoverable : std::exception {
  explicit SurfaceRecoverable(int code) : rc{code} {}
  int rc{0};
  [[nodiscard]] const char *what() const noexcept override {
    return "host surface recoverable";
  }
};

class HostSurfacePotential : public Potential {
public:
  HostSurfacePotential(eon_relax_surface_fn fn, void *user,
                       std::uint64_t epoch)
      : Potential(PotType::UNKNOWN), fn_{fn}, user_{user}, epoch_{epoch} {}

  [[nodiscard]] bool supportsBatchEvaluation() const noexcept override {
    return true;
  }
  [[nodiscard]] unsigned long long surfaceEpoch() const noexcept override {
    return epoch_;
  }
  [[nodiscard]] int lastStatus() const noexcept { return last_; }
  void setEpoch(std::uint64_t epoch) { epoch_ = epoch; }
  void bindPath(const std::vector<std::shared_ptr<Matter>> &path) {
    path_ = path;
  }

  void force(long nAtoms, const double *positions, const int *atomicNrs,
             double *forces, double *energy, double *variance,
             const double *box) override {
    const double *posv[] = {positions};
    const int *nrv[] = {atomicNrs};
    double *frcv[] = {forces};
    const double *boxv[] = {box};
    double e{0};
    double var{0};
    count_calls_ = false;
    forceBatch(1, nAtoms, posv, nrv, frcv, &e, &var, boxv);
    count_calls_ = true;
    if (energy) {
      *energy = e;
    }
    if (variance) {
      *variance = var;
    }
  }

  void forceBatch(long nSystems, long nAtoms, const double *const *positions,
                  const int *const *atomicNrs, double *const *forces,
                  double *energies, double *variances,
                  const double *const *boxes) override {
    if (!fn_) {
      throw std::runtime_error("HostSurfacePotential: NULL surface");
    }
    const long npos = 3 * nAtoms * nSystems;
    const long nbox = 9 * nSystems;
    std::vector<double> pos(static_cast<size_t>(npos), 0.0);
    std::vector<double> frc(static_cast<size_t>(npos), 0.0);
    std::vector<double> box(static_cast<size_t>(nbox), 0.0);
    std::vector<int32_t> z(static_cast<size_t>(nAtoms), 0);
    std::vector<int64_t> ids(static_cast<size_t>(nSystems), 0);
    if (atomicNrs && atomicNrs[0]) {
      for (long a = 0; a < nAtoms; ++a) {
        z[static_cast<size_t>(a)] = atomicNrs[0][a];
      }
    }
    for (long s = 0; s < nSystems; ++s) {
      int64_t id = EON_RELAX_IMAGE_NONE;
      if (positions && positions[s]) {
        for (size_t i = 0; i < path_.size(); ++i) {
          if (path_[i] && path_[i]->getPositions().data() == positions[s]) {
            id = static_cast<int64_t>(i);
            break;
          }
        }
      }
      ids[static_cast<size_t>(s)] = id;
      if (positions && positions[s]) {
        for (long k = 0; k < 3 * nAtoms; ++k) {
          pos[static_cast<size_t>(s * 3 * nAtoms + k)] = positions[s][k];
        }
      }
      if (boxes && boxes[s]) {
        for (long k = 0; k < 9; ++k) {
          box[static_cast<size_t>(s * 9 + k)] = boxes[s][k];
        }
      }
    }
    std::vector<double> ebuf(static_cast<size_t>(nSystems), 0.0);
    std::vector<double> vbuf(static_cast<size_t>(nSystems), 0.0);
    std::uint64_t new_epoch = epoch_;
    eon_relax_surface_request_t req{};
    req.version.major = EON_RELAX_ABI_MAJOR;
    req.version.minor = EON_RELAX_ABI_MINOR;
    req.n_images = static_cast<int64_t>(nSystems);
    req.n_atoms = static_cast<int64_t>(nAtoms);
    req.positions = pos.data();
    req.atomic_nrs = z.data();
    req.boxes = box.data();
    req.image_ids = ids.data();
    req.energies = ebuf.data();
    req.forces = frc.data();
    req.variances = vbuf.data();
    req.epoch_out = &new_epoch;
    last_ = fn_(user_, &req);
    if (last_ < 0) {
      throw std::runtime_error("HostSurfacePotential: surface failed");
    }
    if (last_ > 0) {
      throw SurfaceRecoverable(last_);
    }
    if (new_epoch != epoch_) {
      epoch_ = new_epoch;
    }
    if (count_calls_) {
      for (long s = 0; s < nSystems; ++s) {
        forceCallCounter++;
        PotRegistry::get().on_force_call(ptype);
      }
    }
    for (long s = 0; s < nSystems; ++s) {
      if (energies) {
        energies[s] = ebuf[static_cast<size_t>(s)];
      }
      if (variances) {
        variances[s] = vbuf[static_cast<size_t>(s)];
      }
      if (forces && forces[s]) {
        for (long k = 0; k < 3 * nAtoms; ++k) {
          forces[s][k] = frc[static_cast<size_t>(s * 3 * nAtoms + k)];
        }
      }
    }
  }

private:
  eon_relax_surface_fn fn_;
  void *user_;
  std::uint64_t epoch_{0};
  int last_{0};
  bool count_calls_{true};
  std::vector<std::shared_ptr<Matter>> path_;
};

} // namespace eonc
