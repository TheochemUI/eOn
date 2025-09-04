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

class ExtPot : public Potential {

public:
  ExtPot(std::shared_ptr<Parameters> p)
      : Potential(p),
        eon_extpot_path{p->extPotPath.c_str()} {};
  ~ExtPot();
  void cleanMemory(void);
  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box) override;

private:
  void passToSystem(long N, const double *R, const int *atomicNrs,
                    const double *box);
  void recieveFromSystem(long N, double *F, double *U);
  const char *eon_extpot_path;
};
