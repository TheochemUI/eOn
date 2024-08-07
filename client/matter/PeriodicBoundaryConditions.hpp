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
#include "client/Eigen.h"
namespace eonc::mat {

/**
 * @brief Base class for Periodic Boundary Condition (PBC).
 *
 * This class defines the interface for applying periodic boundary conditions.
 * Derived classes must implement the `operator()` method to handle specific
 * PBC algorithms.
 */
class PeriodicBoundaryCondition {
public:
  /**
   * @brief Construct a new Periodic Boundary Condition object.
   *
   * @param cell The simulation cell matrix.
   * @param cellInverse The inverse of the simulation cell matrix.
   */
  PeriodicBoundaryCondition(const Matrix3S &cell)
      : cell(cell),
        cellInverse(cell.inverse()) {}

  virtual ~PeriodicBoundaryCondition() = default;
  virtual AtomMatrix operator()(const AtomMatrix &) const = 0;

protected:
  const Matrix3S cell;
  const Matrix3S cellInverse;
};

/**
 * @brief PBC Algorithm A1 Implementation.
 *
 * Uses an implicit double-to-integer cast operation for efficient remainder
 * calculation.
 *
 * @note Wrapped Range :: [0, 1).
 *
 * Reference: Dieters, Algorithm A1.
 */
class PBC_A1 : public PeriodicBoundaryCondition {
public:
  using PeriodicBoundaryCondition::PeriodicBoundaryCondition;
  AtomMatrix operator()(const AtomMatrix &) const override;
};

/**
 * @brief PBC Algorithm B1 Implementation.
 *
 * Extends A1 by performing the "cast" step twice for a wider working range.
 *
 * @note Wrapped Range :: [0, 1).
 *
 * Reference: Dieters, Algorithm B1.
 */
class PBC_B1 : public PeriodicBoundaryCondition {
public:
  using PeriodicBoundaryCondition::PeriodicBoundaryCondition;
  AtomMatrix operator()(const AtomMatrix &) const override;
};

/**
 * @brief PBC Algorithm B2 Implementation.
 *
 * Uses a conditional clause to adjust coordinates within the primary simulation
 * box.
 *
 * @note Wrapped Range :: [-0.5, 0.5).
 *
 * Reference: Dieters, Algorithm B2.
 */
class PBC_B2 : public PeriodicBoundaryCondition {
public:
  using PeriodicBoundaryCondition::PeriodicBoundaryCondition;
  AtomMatrix operator()(const AtomMatrix &) const override;
};

/**
 * @brief PBC Algorithm B3 Implementation.
 *
 * Simplifies B2 by disregarding the sign of the coordinate difference.
 *
 * @note Wrapped Range :: [-0.5, 0.5).
 *
 * Reference: Dieters, Algorithm B3.
 */
class PBC_B3 : public PeriodicBoundaryCondition {
public:
  using PeriodicBoundaryCondition::PeriodicBoundaryCondition;
  AtomMatrix operator()(const AtomMatrix &) const override;
};

/**
 * @brief PBC Algorithm B4 Implementation.
 *
 * Sidesteps negative distances by adding a multiple of the box length.
 *
 * @note Wrapped Range :: [-3.5b, +∞).
 *
 * Reference: Dieters, Algorithm B4.
 */
class PBC_B4 : public PeriodicBoundaryCondition {
public:
  using PeriodicBoundaryCondition::PeriodicBoundaryCondition;
  AtomMatrix operator()(const AtomMatrix &) const override;
};

/**
 * @brief PBC Algorithm C1 Implementation.
 *
 * Uses the `remainder` function for periodic boundary condition adjustment.
 *
 * @note Wrapped Range :: [-0.5, 0.5).
 *
 * Reference: Dieters, Algorithm C1.
 */
class PBC_C1 : public PeriodicBoundaryCondition {
public:
  using PeriodicBoundaryCondition::PeriodicBoundaryCondition;
  AtomMatrix operator()(const AtomMatrix &) const override;
};

/**
 * @brief PBC Algorithm C5 Implementation.
 *
 * Uses a fast double-to-integer cast for periodic boundary condition
 * adjustment.
 *
 * @note Wrapped Range :: [-0.5, 0.5).
 *
 * Reference: Dieters, Algorithm C5.
 */
class PBC_C5 : public PeriodicBoundaryCondition {
public:
  using PeriodicBoundaryCondition::PeriodicBoundaryCondition;
  AtomMatrix operator()(const AtomMatrix &) const override;
};

/**
 * @brief PBC Algorithm C6 Implementation.
 *
 * Simplifies C5 by disregarding the sign of the coordinate difference.
 *
 * @note Wrapped Range :: [-0.5, 0.5).
 *
 * Reference: Dieters, Algorithm C6.
 */
class PBC_C6 : public PeriodicBoundaryCondition {
public:
  using PeriodicBoundaryCondition::PeriodicBoundaryCondition;
  AtomMatrix operator()(const AtomMatrix &) const override;
};

/**
 * @brief PBC eOn SVN Implementation.
 *
 * Older SVN variant extracted from Matter
 *
 * @note Wrapped Range :: [-0.5, 0.5).
 *
 * Reference: Dieters, Algorithm C6.
 */
class PBC_ESVN : public PeriodicBoundaryCondition {
public:
  using PeriodicBoundaryCondition::PeriodicBoundaryCondition;
  AtomMatrix operator()(const AtomMatrix &) const override;
};

} // namespace eonc::mat
