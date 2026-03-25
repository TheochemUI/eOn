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
#include "Matter.h"
#include "Parameters.h"

#include <vector>

#include "Eigen.h"
#include "EonLogger.h"

namespace eonc {

/** Functionality relying on the conjugate gradients algorithm. The object is
 * capable of minimizing an Matter object or modified forces being passed in.*/
class BondBoost {

public:
  /** Constructor to be used when a structure is minimized.
  @param[in]   *matter        Pointer to the Matter object to be relaxed.
  @param[in]   parameters     Reference to the Parameter object containing the
  runtime parameters.*/
  BondBoost(Matter *matt, const Parameters &params);
  ~BondBoost(); ///< Destructor.

  void initialize();
  double boost();

private:
  Matrix<double, Eigen::Dynamic, 1> Rmdsteps();
  long BondSelect();
  double Booststeps();
  long nAtoms{0}; ///< Number of free coordinates.
  Matter *matter{
      nullptr}; ///< Pointer to atom object \b outside the scope of the class.
  const Parameters &parameters; ///< Reference to a structure outside the scope
                                ///< of the class containing runtime parameters.
  std::vector<long> BAList;
  std::vector<long> RAList;
  std::vector<long> TABAList;
  std::vector<long> BBAList;
  std::vector<double> Epsr_Q;
  Matrix<double, Eigen::Dynamic, 1>
      TABLList; // EquilibriumTaggedAtomInvolvedBondLengthList;
  Matrix<double, Eigen::Dynamic, 1> EBBLList; // EquilibriumBoostBondLengthList
  Matrix<double, Eigen::Dynamic, 1> CBBLList; // CurrentBoostBondLengthList
  long nBAs{0};
  long nRAs{0};
  long nTABs{0};
  long nReg{0};
  long nBBs{0};
  eonc::log::Scoped log;
};

class Hyperdynamics {
public:
  //            static const string NONE;
  //            static const string BOND_BOOST;
  static const char NONE[];
  static const char BOND_BOOST[];
};

} // namespace eonc

using eonc::BondBoost;
using eonc::Hyperdynamics;
