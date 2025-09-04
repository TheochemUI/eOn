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

/** @file
      @brief Morse potential for platinum
      @author Anonymous (possibly A. Pedersen or G. Henkelman), revision: Jean
   Claude C. Berthet
      @date Unknown, revision: 2010, University of Iceland
      */
// #include "LJBinary.h"
#include "../../Potential.h"
#include <cmath>

class Morse : public Potential {
public:
  Morse(std::shared_ptr<Parameters> params)
      : Potential(params),
        De_{0.7102},
        a_{1.6047},
        re_{2.8970},
        cutoff_{9.5} {
    setParameters(De_, a_, re_, cutoff_);
  };
  // Parameters De in eV, a in Angstroms, re in Angstroms, cutoff in Angstroms
  // Morse(double re, double De, double a, double cutoff);
  void force(long N, const double *R, const int *, double *F, double *U,
             double *variance, const double *box) override;
  void setParameters(double De, double a, double re, double cutoff);

private:
  void morse(double r, double &energy, double &force);
  double De_;
  double a_;
  double re_;
  double cutoff_;
  double energyCutoff_;
};
