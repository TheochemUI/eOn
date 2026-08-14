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

#include <filesystem>
#include <string>

class ExtPot : public Potential {

public:
  ExtPot(const Parameters &p);
  ~ExtPot();
  void cleanMemory(void);
  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box) override;

private:
  void passToSystem(long N, const double *R, const int *atomicNrs,
                    const double *box);
  void recieveFromSystem(long N, double *F, double *U);
  //!< Create the directory the external command runs in and clear any
  //!< exchange file left in it.
  void prepareExchangeDir();
  // Owned by value: the Parameters this was built from may go out of scope
  // long before the potential does.
  std::string eon_extpot_path;
  //!< The directory this potential runs the external command in. The
  //!< exchange files keep their documented names and live here, so two
  //!< clients started in one directory cannot read each other's numbers.
  std::filesystem::path exchangeDir;
  //!< Set while a call is in flight, so a call that fails leaves its
  //!< exchange files in place to be looked at.
  bool retainExchangeDir{false};
};
