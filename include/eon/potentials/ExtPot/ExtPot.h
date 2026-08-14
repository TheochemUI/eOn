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

#include "eon/Potential.h"

#include <string>

class ExtPot : public Potential {

public:
  ExtPot(const Parameters &p)
      : Potential(p),
        eon_extpot_path{p.potential_options.extPotPath} {};
  ~ExtPot();
  void cleanMemory(void);
  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box) override;

private:
  void passToSystem(long N, const double *R, const int *atomicNrs,
                    const double *box);
  void recieveFromSystem(long N, double *F, double *U);
  // Owned by value: the Parameters this was built from may go out of scope
  // long before the potential does.
  std::string eon_extpot_path;
};
