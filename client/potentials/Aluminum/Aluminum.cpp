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
//-----------------------------------------------------------------------------------
#include "Aluminum.h"
#include <assert.h>

namespace eonc {

void Aluminum::forceImpl(const ForceInput &params, ForceOut *efvd) {
#ifdef EON_CHECKS
  eonc::pot::checkParams(params);
  eonc::pot::zeroForceOut(params.nAtoms, efvd);
#endif
  const long int N = params.nAtoms;
  force_(&N, params.pos, efvd->F, &efvd->energy, &params.box[0], &params.box[4],
         &params.box[8]);
  return;
}

} // namespace eonc
