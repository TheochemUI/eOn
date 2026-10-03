// The MPI half of RgpotPot's calculator groups: the construction agreement
// before the calculator split and the grouped exit handler. It is a library
// of its own so that librgpot_pot, and with it every eonclient, carries no
// MPI dependency; RGPotEngine loads it with dlopen only when the process was
// started under an MPI launcher. rgpot loads its own MPI side the same way.
#include "eon/potentials/Rgpot/RgpotGroupMpi.h"

#include <mpi.h>

#include <algorithm>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "rgpot/CalculatorGroup.hpp"

namespace {

// Every rank calls this. A rank that failed locally still enters, so nobody
// is left in MPI_Comm_split. The shared message is the lowest failing rank's
// text. Returns 1 when every rank succeeded.
int agree(const char *local, char *shared, std::size_t cap) {
  int inited = 0;
  MPI_Initialized(&inited);
  if (!inited)
    MPI_Init(nullptr, nullptr);
  ::rgpot::finalizeMpiAtExit();
  const std::size_t localLen = local == nullptr ? 0 : std::strlen(local);
  const int ok = localLen == 0 ? 1 : 0;
  int all_ok = 0;
  MPI_Allreduce(&ok, &all_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  if (cap > 0)
    shared[0] = '\0';
  if (all_ok)
    return 1;
  int rank = 0;
  int size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  const int mine = ok ? size : rank;
  int owner = 0;
  MPI_Allreduce(&mine, &owner, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  int len = 0;
  if (rank == owner)
    len = static_cast<int>(std::min(localLen, std::size_t{4095}));
  MPI_Bcast(&len, 1, MPI_INT, owner, MPI_COMM_WORLD);
  std::vector<char> buf(static_cast<std::size_t>(len) + 1, '\0');
  if (rank == owner && len > 0)
    std::memcpy(buf.data(), local, static_cast<std::size_t>(len));
  if (len > 0)
    MPI_Bcast(buf.data(), len, MPI_CHAR, owner, MPI_COMM_WORLD);
  if (cap > 0) {
    const std::size_t n = std::min(static_cast<std::size_t>(len), cap - 1);
    std::memcpy(shared, buf.data(), n);
    shared[n] = '\0';
  }
  return 0;
}

// Runs before rgpot's atexit handler and ends in _Exit, so rgpot's handler
// never runs. It honours rgpot's abort request itself: after a failed
// engine call the peers can sit in a collective, and MPI_Finalize would
// wait for them until the walltime kill.
void hardExit(int status, void *) {
  int inited = 0;
  int finalized = 0;
  MPI_Initialized(&inited);
  MPI_Finalized(&finalized);
  if (inited && !finalized) {
    if (::rgpot::mpiAbortRequested()) {
      std::fflush(nullptr);
      MPI_Abort(MPI_COMM_WORLD, status != 0 ? status : 1);
    }
    MPI_Finalize();
  }
  std::fflush(nullptr);
  std::_Exit(status);
}

const EonRgpotGroupMpi kApi{EON_RGPOT_GROUP_MPI_VERSION, agree, hardExit};

} // namespace

extern "C" const EonRgpotGroupMpi *eon_rgpot_group_mpi_v1() { return &kApi; }
