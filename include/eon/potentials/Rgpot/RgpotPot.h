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

#include "eon/Potential.h"
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
  /// spread over the groups: system j runs on group j mod G, then every
  /// rank receives every result, so all ranks keep the same band.
  [[nodiscard]] bool supportsBatchEvaluation() const noexcept override;
  void forceBatch(long nSystems, long nAtoms, const double *const *positions,
                  const int *const *atomicNrs, double *const *forces,
                  double *energies, double *variances,
                  const double *const *boxes) override;
  /// NWChem molecular SCF does not support PBC; CPMD may be periodic.
  [[nodiscard]] bool requiresIsolatedMoleculeLayout() const noexcept override {
    return backend_ == "nwchemc";
  }
  [[nodiscard]] const std::string &backend() const noexcept { return backend_; }
  [[nodiscard]] bool engineAvailable() const;

private:
  // Under mpirun (an MPI build of rgpot with cpmdc) only world rank 0 runs
  // eOn. It broadcasts each force request; the other ranks serve requests
  // from the constructor until the driver sends a stop, then exit.
  [[noreturn]] void serveWorker();
  void computeSingle(long N, const double *R, const int *atomicNrs, double *F,
                     double *U, const double *box);
  void computeBatch(long nSystems, long nAtoms, const double *const *positions,
                    const int *const *atomicNrs, double *const *forces,
                    double *energies, const double *const *boxes);
  void sendStop();
  void releaseWorkersAtExit();

  std::unique_ptr<RGPotEngine> impl_;
  std::string backend_;
  bool driver_{true};
  bool stopped_{false};
};
