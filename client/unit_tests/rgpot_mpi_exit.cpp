/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
*/
// Three mpirun launches of one binary.
//   params: rank 1 cannot read params_path. Both ranks must throw that
//           error. A rank left in MPI_Comm_split never prints rank=0.
//   fault:  rank 1's engine throws engine-rank1. Rank 0 must print it
//           and exit non-zero. Workers stop with status 0.
//   abort:  rank 0 requests rgpot's abort-at-exit, as a failed engine
//           call does, and exits while rank 1 waits in the worker
//           broadcast. The exit handler must abort the world; an
//           MPI_Finalize there waits for rank 1 until the timeout.
//   single: one structure through force() on two calculators. The
//           engine agrees on errors across the world after every call,
//           so both calculators must call it; a calculator left out
//           holds the other in that agreement until the timeout.
//   uneven: three structures through forceBatchOwned on two
//           calculators, two owned by one and one by the other. Both
//           calculators must make the same number of engine calls.

#include "eon/Parameters.h"
#include "eon/Potential.h"
#include "rgpot/CalculatorGroup.hpp"

#include <cstdlib>
#include <iostream>
#include <span>
#include <string>

namespace {

bool env_nonempty(const char *name) {
  const char *v = std::getenv(name);
  return v != nullptr && v[0] != '\0';
}

int world_rank() {
  const char *r = std::getenv("OMPI_COMM_WORLD_RANK");
  if (r == nullptr || r[0] == '\0')
    r = std::getenv("PMI_RANK");
  if (r == nullptr || r[0] == '\0')
    return -1;
  return std::atoi(r);
}

} // namespace

int main(int argc, char **argv) {
  if (argc != 2)
    return 2;
  const std::string mode = argv[1];
  if (mode != "params" && mode != "fault" && mode != "abort" &&
      mode != "single" && mode != "uneven")
    return 2;
  if (!(env_nonempty("CPMDC_LIBRARY") || env_nonempty("RGPOT_CPMDC_ENGINE") ||
        env_nonempty("RGPOT_CPMD_ENGINE"))) {
    return 77;
  }
  const int rank = world_rank();
  if (rank < 0) {
    std::cerr << "rgpot_mpi_exit: not started under mpirun\n";
    return 77;
  }

  if (mode == "params") {
    if (rank == 1)
      setenv("RGPOT_PARAMS_PATH", "missing-params.bin", 1);
    else
      unsetenv("RGPOT_PARAMS_PATH");
  } else if (mode == "fault" && rank == 1) {
    setenv("RGPOT_FORCE_FAIL", "engine-rank1", 1);
  } else {
    unsetenv("RGPOT_FORCE_FAIL");
  }

  eonc::Parameters params{};
  eonc::ParametersLoadAccess::potential_options(params).potential =
      eonc::PotType::RGPOT;
  eonc::ParametersLoadAccess::rgpot_options(params).backend = "cpmdc";
  eonc::ParametersLoadAccess::rgpot_options(params).functional = "BLYP";
  // single, uneven and fault run real SCFs, so a small cell at a low
  // cutoff keeps each call to seconds.
  const bool scf = mode == "single" || mode == "uneven" || mode == "fault";
  eonc::ParametersLoadAccess::rgpot_options(params).cutoff_ry =
      scf ? 30.0 : 70.0;
  eonc::ParametersLoadAccess::rgpot_options(params).charge = 0;
  eonc::ParametersLoadAccess::rgpot_options(params).multiplicity = 1;
  if (env_nonempty("CPMDC_LIBRARY"))
    eonc::ParametersLoadAccess::rgpot_options(params).engine_path =
        std::getenv("CPMDC_LIBRARY");
  // Two ranks, one rank per calculator, so the fault sits on group 1.
  eonc::ParametersLoadAccess::rgpot_options(params).ranks_per_image =
      mode == "params" ? 0 : 1;

  try {
    auto pot = eonc::helpers::sharePotential(
        eonc::helpers::makePotential(eonc::PotType::RGPOT, params));
    if (mode == "abort") {
      // Only the driver returns from the grouped constructor. std::exit
      // keeps pot alive, so no stop message reaches rank 1.
      ::rgpot::abortMpiAtExit();
      std::cerr << "rank=" << rank << " abort requested\n";
      std::exit(3);
    }
    if (mode == "single" || mode == "uneven") {
      // N2, closed shell, at three bond lengths: the SCF converges in a
      // few iterations in a 6 A cell.
      const long n = mode == "single" ? 1 : 3;
      const double bond[3] = {1.10, 1.12, 1.08};
      double R[3][6] = {};
      int Z[3][2] = {{7, 7}, {7, 7}, {7, 7}};
      double F[3][6] = {};
      double U[3] = {};
      double box[9] = {6, 0, 0, 0, 6, 0, 0, 0, 6};
      for (long j = 0; j < 3; j++) {
        R[j][0] = 2.6;
        R[j][1] = R[j][2] = 3.0;
        R[j][3] = 2.6 + bond[j];
        R[j][4] = R[j][5] = 3.0;
      }
      if (mode == "single") {
        double var = 0.0;
        pot->force(std::span<const double>(R[0], 6),
                   std::span<const int>(Z[0], 2), std::span<double>(F[0], 6), U,
                   &var, std::span<const double>(box, 9));
      } else {
        const double *pos[3] = {R[0], R[1], R[2]};
        const int *nrs[3] = {Z[0], Z[1], Z[2]};
        double *frc[3] = {F[0], F[1], F[2]};
        const double *bx[3] = {box, box, box};
        long owners[3] = {0, 1, 2};
        pot->forceBatchOwned(n, 2, pos, nrs, frc, U, nullptr, bx, owners);
      }
      std::cerr << "rank=" << rank << " " << mode << " done E0=" << U[0]
                << "\n";
      return 0;
    }
    if (mode != "fault") {
      std::cerr << "rank=" << rank << " params constructed\n";
      return 1;
    }
    // N2 in a 6 A cell, as in single and uneven.
    double R[6] = {2.6, 3.0, 3.0, 3.7, 3.0, 3.0};
    int Z[2] = {7, 7};
    double F[6] = {};
    double U = 0.0;
    double box[9] = {6, 0, 0, 0, 6, 0, 0, 0, 6};
    const double *pos = R;
    const int *nrs = Z;
    double *frc = F;
    const double *bx = box;
    long owner = 1;
    try {
      pot->forceBatchOwned(1, 2, &pos, &nrs, &frc, &U, nullptr, &bx, &owner);
    } catch (const std::exception &ex) {
      std::cerr << "rank=" << rank << " fault " << ex.what() << "\n";
      return 1;
    }
    std::cerr << "rank=" << rank << " fault did not throw\n";
    return 1;
  } catch (const std::exception &ex) {
    std::cerr << "rank=" << rank << " params-agreed " << ex.what() << "\n";
    return 1;
  }
}
