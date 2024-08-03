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
#include "../../Potential.h"
#include <filesystem>
#include <future>

namespace eonc {

class VASP : public Potential<VASP> {

public:
  VASP() {
    vaspRunCount++;
    // Deleting leftovers from previous run
    const std::vector<std::string> filesToRemove = {
        "TMPCAR",   "CHG",    "CHGCAR", "CONTCAR", "DOSCAR",
        "EIGENVAL", "IBZKPT", "NEWCAR", "FU",      "OSZICAR",
        "OUTCAR",   "PCDAT",  "POSCAR", "WAVECAR", "XDATCAR"};

    for (const auto &file : filesToRemove) {
      std::filesystem::remove(file);
    }
  }
  ~VASP() { cleanMemory(); }
  void cleanMemory(void);
  void forceImpl(const ForceInput &, ForceOut *) override final;

private:
  void writePOSCAR(long N, const double *R, const size_t *atomicNrs,
                   const double *box);
  void readFU(long N, double *F, double *U);
  void spawnVASP();
  bool vaspRunning();
  std::future<void> vaspProcess;
  static bool firstRun;
  static long vaspRunCount;
  static pid_t vaspPID;
};
} // namespace eonc
