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
#include "EonLogger.h"

#include "FiniteDifference.h"
#include "Matter.h"
#include "Parameters.h"

#include <vector>

namespace eonc {

class Hessian {
public:
  Hessian(const Parameters &params, Matter *matter);
  ~Hessian() = default;

  MatrixXd getHessian(Matter *matterIn, const VectorXi &atomsIn);
  VectorXd getFreqs(Matter *matterIn, const VectorXi &atomsIn);
  //    VectorXd getModes(Matter *matterIn, VectorXi atomsIn);
  VectorXd removeZeroFreqs(const VectorXd &freqs);

private:
  Matter *matter;
  const Parameters &parameters;

  MatrixXd hessian;
  VectorXd freqs;

  VectorXi atoms;
  bool calculate();
  bool finalizeHessian(int size);
  bool calculateColored(double cutoff, double dr, FdScheme scheme);
  bool calculateSerial(double dr, FdScheme scheme);
  eonc::log::Scoped log;
};

/// Greedy coloring of an undirected graph. ``adj[v]`` lists neighbors of
/// ``v`` (self-loops ignored). Colors are dense from 0 in vertex order.
std::vector<int>
greedyColorCutoffGraph(const std::vector<std::vector<int>> &adj);

/// Conflict graph on ``atoms``: two mobile atoms share an edge when their
/// closed cutoff neighborhoods intersect, then greedy-colored. Empty when
/// the neighbor list cannot be built.
std::vector<int> colorMobileCutoffGraph(const Matter &matter,
                                        const VectorXi &atoms, double cutoff);

} // namespace eonc
