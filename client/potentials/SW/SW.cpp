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

#include "SW.h"
namespace eonc {
void SW::forceImpl(const ForceInput &fip, ForceOut *efvd) {
#ifdef EON_CHECKS
  eonc::pot::checkParams(fip);
  eonc::pot::zeroForceOut(fip.nAtoms, efvd);
#endif
  const long int N = fip.nAtoms;
  sw_(&N, fip.pos, efvd->F, &efvd->energy, &fip.box[0], &fip.box[4],
      &fip.box[8]);
  return;
}

} // namespace eonc
