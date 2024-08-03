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

// serves as an interface between emt potentials provided by CamposASE and
// dynamics provided by EON

#pragma once
#include "Asap/Atoms.h"
#include "Asap/EMT.h"
#include "Asap/EMTRasmussenParameterProvider.h"
#include "Asap/SuperCell.h"

#include "../../Potential.h"
namespace eonc {
/** EMT potential. Inspect the EMT_parms.h to see what the EMT potential is
 * hardcoded to describe.*/
class EffectiveMediumTheory : public Potential<EffectiveMediumTheory> {
private:
  size_t numberOfAtoms;
  bool periodicity[3];
  std::unique_ptr<SuperCell> SuperCellObj;
  std::unique_ptr<Atoms> AtomsObj;
  std::unique_ptr<EMTRasmussenParameterProvider> EMTRasmussenParams;
  std::unique_ptr<EMT> EMTObj;
  bool useEMTRasmussen, usePBC;
  void initialize(const ForceInput &fip);

public:
  EffectiveMediumTheory(bool useEMTRasmussen, bool usePBC)
      : numberOfAtoms{0},
        useEMTRasmussen{useEMTRasmussen},
        usePBC{usePBC} {
    if (usePBC == false) {
      throw std::invalid_argument(
          "EMT should have periodic boundary conditions in all directions");
    }
    periodicity[0] = true;
    periodicity[1] = true;
    periodicity[2] = true;
  }
  void forceImpl(const ForceInput &, ForceOut *) override final;
};

} // namespace eonc
