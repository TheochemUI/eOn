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

class IIRAResource;

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

  /// Production: process-default IRAResource::instance().
  IRACompare();
  /// Test seam: injected resource, no process-default load.
  IRACompare(IIRAResource &resource);

  /// Match two structures using CShDA + SVD (optimal rotation + assignment).
  /// Both structures must have the same atom types.
  MatchResult match(const Matter &m1, const Matter &m2, double distThreshold);

  /// Same as match(), from packed (n,3) row-major coordinates and Z arrays.
  /// Used by the Python server so rot_match does not need a dummy Potential.
  MatchResult matchArrays(int nat1, const int *typ1, const double *pos1,
                          int nat2, const int *typ2, const double *pos2,
                          double distThreshold);

  /// Atom assignment under periodic boundary conditions (CShDA only, no
  /// rotation/SVD). Returns permutation and per-atom distances.
  MatchResult matchPBC(const Matter &m1, const Matter &m2,
                       double distThreshold);

  /// Find all symmetry operations of a structure (SOFI algorithm).
  SymmetryResult findSymmetry(const Matter &m, double threshold,
                              bool prescreenIh = true);

  /// Rigid-align + permute reactant onto product. Returns the match
  /// (error != 0 means reactant was left unchanged).
  MatchResult alignReactantToProduct(Matter &reactant, const Matter &product,
                                     double distThreshold);

private:
  IIRAResource &resource_;
};

} // namespace eonc
