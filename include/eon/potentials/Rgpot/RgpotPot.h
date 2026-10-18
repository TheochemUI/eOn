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

#include "eon/CalculatorGroupUse.h"
#include "eon/Potential.h"
#include "eon/potentials/Rgpot/GroupSchedule.h"
#include <cstdint>
#include <memory>
#include <string>

class RGPotEngine;

/**
 * Potential backed by rgpot NWChemPot / CPMDPot (in-process dlopen of
 * libnwchemc / libcpmdc). Not potserv RPC. Configure via [RgpotPot] INI;
 * energy eV, forces eV/Angstrom.
 */
class RgpotPot final : public eonc::Potential {
public:
  explicit RgpotPot(const eonc::Parameters &p);
  ~RgpotPot() override;

  RgpotPot(const RgpotPot &) = delete;
  RgpotPot &operator=(const RgpotPot &) = delete;

  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box) override;

  [[nodiscard]] bool isThreadSafe() const noexcept override { return false; }

  /// With cpmdc calculator groups ([RgpotPot] ranks_per_image), a batch is
  /// spread over the groups by eonc::GroupSchedule: a system stays on the
  /// group that holds its orbitals, at most ceil(M / G) systems per group,
  /// then every rank receives every result, so all ranks keep the same
  /// band.
  [[nodiscard]] bool supportsBatchEvaluation() const noexcept override;
  void forceBatch(long nSystems, long nAtoms, const double *const *positions,
                  const int *const *atomicNrs, double *const *forces,
                  double *energies, double *variances,
                  const double *const *boxes) override;
  /// owners[j] is system j's identity (a NEB image, a ring bead): it names
  /// the orbitals the engine starts from and keeps the system on one group
  /// across steps. Without it, or below zero, the batch position stands in.
  void forceBatchOwned(long nSystems, long nAtoms,
                       const double *const *positions,
                       const int *const *atomicNrs, double *const *forces,
                       double *energies, double *variances,
                       const double *const *boxes, const long *owners) override;
  /// NWChem molecular SCF does not support PBC; CPMD may be periodic.
  [[nodiscard]] bool requiresIsolatedMoleculeLayout() const noexcept override {
    return backend_ == "nwchemc";
  }
  [[nodiscard]] const std::string &backend() const noexcept { return backend_; }
  [[nodiscard]] bool engineAvailable() const;
  /// Calculator-group accounting on the driver. The per-group columns are
  /// filled when the driver stops the workers; the driver prints them then.
  [[nodiscard]] const eonc::CalculatorGroupUse &groupUse() const noexcept {
    return use_;
  }

private:
  // Under mpirun (an MPI build of rgpot with cpmdc) only world rank 0 runs
  // eOn. It broadcasts each force request; the other ranks serve requests
  // from the constructor until the driver sends a stop, then exit.
  [[noreturn]] void serveWorker();
  void computeSingle(long N, const double *R, const int *atomicNrs, double *F,
                     double *U, const double *box, int group);
  void computeBatch(long nSystems, long nAtoms, const double *const *positions,
                    const int *const *atomicNrs, double *const *forces,
                    double *energies, const double *const *boxes,
                    const std::int64_t *owners, const std::int64_t *route);
  void sendStop();
  void stopAndDrop();
  void releaseWorkersAtExit();
  // Collective over every rank: each group's first rank shares its busy
  // seconds and system count.
  void exchangeGroupUse();
  bool evaluate(long N, const double *R, const int *atomicNrs, double *F,
                double *U, const double *box, std::string &error);

  std::unique_ptr<RGPotEngine> impl_;
  std::string backend_;
  bool driver_{true};
  bool stopped_{false};
  bool dropped_{false};
  bool acked_{false};
  // This rank's seconds inside engine calls and the systems it evaluated.
  double busy_{0.0};
  double systemsDone_{0.0};
  eonc::CalculatorGroupUse use_;
  // Driver only: which group evaluates which system.
  std::unique_ptr<eonc::GroupSchedule> schedule_;
};
