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

#include "eon/potentials/MPIPot/MPIPot.h"
#include "eon/Parameters.h"

#include <array>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <mpi.h>
#include <stdexcept>
#include <string>
#include <thread>

namespace {

void mpi_check(int rc, const char *what) {
  if (rc != MPI_SUCCESS) {
    throw std::runtime_error(std::string(what) + " failed");
  }
}

} // namespace

MPIPot::MPIPot(const eonc::Parameters &p)
    : eonc::Potential(p) {
  potentialRank = p.potential_options().MPIPotentialRank;
  poll_period = p.potential_options().MPIPollPeriod;
}

void MPIPot::cleanMemory() {}

MPIPot::~MPIPot() { cleanMemory(); }

void MPIPot::force(long N, const double *R, const int *atomicNrs, double *F,
                   double *U, double *variance, const double *box) {
  if (variance != nullptr) {
    *variance = 0.0;
  }
  // The peer reads 1024 MPI_INT values from this long buffer. The count and
  // datatype stay MPI_INT so the bytes on the wire do not change.
  std::array<long, 1024> icwd{};
  const std::string cwd = std::filesystem::current_path().string();
  if (cwd.size() >= icwd.size()) {
    throw std::runtime_error(
        "working directory path exceeds the MPI potential buffer");
  }
  for (std::size_t i = 0; i < cwd.size(); ++i) {
    icwd[i] = static_cast<long>(static_cast<unsigned char>(cwd[i]));
  }
  int pbc = 1;
  int failed = 0;
  const int intn = static_cast<int>(N);
  mpi_check(MPI_Send(&intn, 1, MPI_INT, potentialRank, 0, MPI_COMM_WORLD),
            "MPI_Send atom count");
  mpi_check(
      MPI_Send(atomicNrs, intn, MPI_INT, potentialRank, 0, MPI_COMM_WORLD),
      "MPI_Send atomic numbers");
  mpi_check(MPI_Send(R, 3 * intn, MPI_DOUBLE, potentialRank, 0, MPI_COMM_WORLD),
            "MPI_Send positions");
  mpi_check(MPI_Send(box, 9, MPI_DOUBLE, potentialRank, 0, MPI_COMM_WORLD),
            "MPI_Send cell");
  mpi_check(MPI_Send(&pbc, 1, MPI_INT, potentialRank, 0, MPI_COMM_WORLD),
            "MPI_Send pbc");
  mpi_check(MPI_Send(icwd.data(), static_cast<int>(icwd.size()), MPI_INT,
                     potentialRank, 0, MPI_COMM_WORLD),
            "MPI_Send working directory");

  if (poll_period > 0.0) {
    int eon_flag = 0;
    // poll_period is divided by 1e6 before the sleep, and the result is
    // truncated to whole microseconds, matching the previous usleep call.
    const auto usec = static_cast<std::uint64_t>(poll_period / 1000000.0);
    mpi_check(MPI_Iprobe(potentialRank, 0, MPI_COMM_WORLD, &eon_flag,
                         MPI_STATUS_IGNORE),
              "MPI_Iprobe");
    while (!eon_flag) {
      if (usec > 0) {
        std::this_thread::sleep_for(std::chrono::microseconds(usec));
      }
      mpi_check(MPI_Iprobe(potentialRank, 0, MPI_COMM_WORLD, &eon_flag,
                           MPI_STATUS_IGNORE),
                "MPI_Iprobe");
    }
  }

  mpi_check(MPI_Recv(&failed, 1, MPI_INT, potentialRank, 0, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE),
            "MPI_Recv status");
  if (failed == 1) {
    throw std::runtime_error("MPI potential reported a failed evaluation");
  }

  mpi_check(MPI_Recv(U, 1, MPI_DOUBLE, potentialRank, 0, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE),
            "MPI_Recv energy");
  mpi_check(MPI_Recv(F, 3 * intn, MPI_DOUBLE, potentialRank, 0, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE),
            "MPI_Recv forces");
}
