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

// serves as an interface between emt potentials provided by CamposASE and
// dynamics provided by eOn

#include "Asap/Atoms.h"
#include "Asap/EMT.h"
#include "Asap/EMTDefaultParameterProvider.h"
#include "Asap/EMTRasmussenParameterProvider.h"
#include "Asap/SuperCell.h"
#include "Asap/Vec.h"

#include "../../Parameters.h"
#include "../../Potential.h"

/** EMT potential. Inspect the EMT_parms.h to see what the EMT potential is
 * hardcoded to describe.*/
class EffectiveMediumTheory : public Potential {

private:
  bool emtRasmussen{false};
  long numberOfAtoms{0};
  bool periodicity[3]{true, true, true};
  Atoms *AtomsObj{nullptr};
  EMTDefaultParameterProvider *EMTParameterObj{nullptr};
  EMT *EMTObj{nullptr};
  SuperCell *SuperCellObj{nullptr};

public:
  EffectiveMediumTheory(const Parameters &p)
      : Potential(p),
        emtRasmussen{p.potential_options.EMTRasmussen} {}
  ~EffectiveMediumTheory() { cleanMemory(); }
  void cleanMemory();

  // To satify interface
  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box);
};
