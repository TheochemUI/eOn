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

/// Optional MPI view of Parameters::potential_options.MPIClientComm.
/// Include this only from TUs that already include mpi.h.

#include "eon/Parameters.h"
#include <cstring>

#ifdef EONMPI
#include "mpi.h"

namespace eonc {

inline void setMpiClientComm(Parameters &p, MPI_Comm comm) {
  static_assert(sizeof(MPI_Comm) <= sizeof(std::uintptr_t),
                "MPI_Comm does not fit in uintptr_t");
  std::uintptr_t raw = 0;
  std::memcpy(&raw, &comm, sizeof(comm));
  p.potential_options.MPIClientComm = raw;
}

inline MPI_Comm getMpiClientComm(const Parameters &p) {
  MPI_Comm comm{};
  std::memcpy(&comm, &p.potential_options.MPIClientComm, sizeof(comm));
  return comm;
}

} // namespace eonc
#endif
