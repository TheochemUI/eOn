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
// serves as an interface between LAMMPS potentials maintained by SANDIA

#pragma once

#include "../../Parameters.h"
#include "../../Potential.h"

#include <filesystem>

/// Thread safety: NOT safe to share an instance across threads.
/// Internally calls std::filesystem::current_path() (process-wide
/// state) and uses a `shell cd` LAMMPS command to pin liblammps's
/// working directory; both are global side effects. liblammps itself
/// also keeps process-global state per LAMMPSObj. Use one LAMMPSPot
/// instance per thread / per NEB image.
class LAMMPSPot : public Potential {

public:
  LAMMPSPot(const Parameters &p);
  ~LAMMPSPot() override;
  void cleanMemory();
  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box) override;

private:
  int lammpsThr{0};
#ifdef EONMPI
  MPI_Comm mpiComm;
#endif
  long numberOfAtoms{0};
  double oldBox[9]{};
  void *LAMMPSObj{nullptr};
  void makeNewLAMMPS(long N, const double *R, const int *atomicNrs,
                     const double *box);
  /// Working directory liblammps reads every relative path from
  /// (in.lammps itself + everything in.lammps references). makeNewLAMMPS
  /// issues `shell cd m_lammps_workdir` to liblammps so the eonclient
  /// CWD doesn't matter.
  ///
  ///   bundle mode (LAMMPSBundlePath set) -> scratch dir under
  ///       temp_directory_path() that LAMMPSBundle::extract() created;
  ///       owned by this instance and removed in the destructor.
  ///   legacy mode (LAMMPSBundlePath empty) -> eonclient's CWD;
  ///       not owned, never removed.
  std::filesystem::path m_lammps_workdir;
  /// Lifetime ownership of m_lammps_workdir. True only in bundle mode.
  bool m_owns_workdir{false};
  bool realunits{false};
};
