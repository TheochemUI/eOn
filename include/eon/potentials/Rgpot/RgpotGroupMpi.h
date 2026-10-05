#pragma once

#include <cstddef>

/// Function table of librgpot_pot_mpi, the MPI half of RgpotPot's
/// calculator groups. librgpot_pot loads it with dlopen and reads it
/// through eon_rgpot_group_mpi_v1, so librgpot_pot itself never links MPI.
#define EON_RGPOT_GROUP_MPI_VERSION 1

extern "C" {

struct EonRgpotGroupMpi {
  int version;
  /// Collective on MPI_COMM_WORLD; initialises MPI when needed. `local`
  /// is this rank's construction error, empty when it succeeded. Returns 1
  /// when every rank succeeded; otherwise `shared` (cap bytes) holds the
  /// lowest failing rank's text on every rank.
  int (*agree)(const char *local, char *shared, std::size_t cap);
  /// on_exit handler: MPI_Abort when rgpot asked for it, else
  /// MPI_Finalize, then _Exit.
  void (*hard_exit)(int status, void *arg);
};

const EonRgpotGroupMpi *eon_rgpot_group_mpi_v1();
}
