/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
*/
// Two mpirun launches of one binary.
//   params: rank 1 cannot read params_path. Both ranks must throw that
//           error. A rank left in MPI_Comm_split never prints rank=0.
//   fault:  rank 1's engine throws engine-rank1. Rank 0 must print it
//           and exit non-zero. Workers stop with status 0.

#include "eon/Parameters.h"
#include "eon/Potential.h"

#include <cstdlib>
#include <iostream>
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
  if (mode != "params" && mode != "fault")
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
  } else if (rank == 1) {
    setenv("RGPOT_FORCE_FAIL", "engine-rank1", 1);
  } else {
    unsetenv("RGPOT_FORCE_FAIL");
  }

  eonc::Parameters params{};
  eonc::ParametersLoadAccess::potential_options(params).potential =
      eonc::PotType::RGPOT;
  eonc::ParametersLoadAccess::rgpot_options(params).backend = "cpmdc";
  eonc::ParametersLoadAccess::rgpot_options(params).functional = "BLYP";
  eonc::ParametersLoadAccess::rgpot_options(params).cutoff_ry = 70.0;
  eonc::ParametersLoadAccess::rgpot_options(params).charge = 0;
  eonc::ParametersLoadAccess::rgpot_options(params).multiplicity = 1;
  if (env_nonempty("CPMDC_LIBRARY"))
    eonc::ParametersLoadAccess::rgpot_options(params).engine_path =
        std::getenv("CPMDC_LIBRARY");
  // Two ranks, one rank per calculator, so the fault sits on group 1.
  eonc::ParametersLoadAccess::rgpot_options(params).ranks_per_image =
      mode == "fault" ? 1 : 0;

  try {
    auto pot = eonc::helpers::makePotential(eonc::PotType::RGPOT, params);
    if (mode != "fault") {
      std::cerr << "rank=" << rank << " params constructed\n";
      return 1;
    }
    double R[3] = {0.0, 0.0, 0.0};
    int Z[1] = {14};
    double F[3] = {};
    double U = 0.0;
    double box[9] = {20, 0, 0, 0, 20, 0, 0, 0, 20};
    const double *pos = R;
    const int *nrs = Z;
    double *frc = F;
    const double *bx = box;
    long owner = 1;
    try {
      pot->forceBatchOwned(1, 1, &pos, &nrs, &frc, &U, nullptr, &bx, &owner);
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
