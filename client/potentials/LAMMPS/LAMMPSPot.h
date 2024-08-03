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
#include "../../Potential.h"

#ifndef EONMPI
#include "fakempi.h"
#endif

namespace eonc {

class LAMMPSPot : public Potential<LAMMPSPot> {

public:
  struct Params final {
    MPI_Comm MPIClientComm{-1};
    int LAMMPSThreads{0};
    std::string suffix{"omp"};
  };

private:
  MPI_Comm _mpicomm;
  int _lmpThreads;
  std::string _suffix;

public:
  LAMMPSPot(const LAMMPSPot::Params &);
  ~LAMMPSPot(void);
  void initialize() {};
  void cleanMemory(void);
  void forceImpl(const ForceInput &, ForceOut *) override final;

private:
  long numberOfAtoms;
  double oldBox[9];
  void *LAMMPSObj;
  void makeNewLAMMPS(long N, const double *R, const size_t *atomicNrs,
                     const double *box);
  bool realunits;
};

} // namespace eonc
