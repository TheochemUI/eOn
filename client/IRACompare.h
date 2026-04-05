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

#include "Eigen.h"
#include "Matter.h"
#include <string>
#include <vector>

namespace eonc {

/// C++ wrapper for the IRA (Iterative Rotations and Assignments) library.
/// Provides structure comparison via CShDA and symmetry detection via SOFI.
class IRACompare {
public:
  struct MatchResult {
    std::vector<int> permutation; // Atom permutation mapping 1->2
    Eigen::Matrix3d rotation;     // Optimal rotation matrix
    Eigen::Vector3d translation;  // Optimal translation vector
    double hausdorffDistance;     // Hausdorff distance after alignment
    int error{0};                 // 0 = success
  };

  struct SymmetryResult {
    std::string pointGroup;
    int nOperations{0};
    std::vector<Eigen::Matrix3d> operations;
    std::vector<std::string> operationLabels;
    std::vector<double> angles;
    std::vector<Eigen::Vector3d> axes;
    int error{0};
  };

  /// Match two structures using CShDA + SVD (optimal rotation + assignment).
  /// Both structures must have the same atom types.
  static MatchResult match(const Matter &m1, const Matter &m2,
                           double distThreshold);

  /// Match with periodic boundary conditions.
  static MatchResult matchPBC(const Matter &m1, const Matter &m2,
                              double distThreshold);

  /// Find all symmetry operations of a structure (SOFI algorithm).
  static SymmetryResult findSymmetry(const Matter &m, double threshold,
                                     bool prescreenIh = true);
};

} // namespace eonc

using eonc::IRACompare;
