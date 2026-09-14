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

namespace eonc {

/// Free (unfixed) atom indices in ascending order.
VectorXi freeAtomIndices(const Matter *matter);

/// Whether atomList means "every free atom" (empty, "all", case-insensitive).
bool atomListMeansAll(const std::string &atomList);

/// PHVA-class mobile set for FD Hessian and matrix-free Krylov (Lanczos /
/// Davidson).
///
/// free/fixed on Matter is the optimizer mask. The mobile set is a
/// (possibly proper) subset of free atoms: those that are *displaced* in
/// finite-difference Hessian-vector products (Li & Jensen, Theor. Chem.
/// Acc. 107, 211 (2002)).
///
/// atomList empty or "all" -> all free atoms (backward-compatible default).
/// Otherwise comma/space-separated 0-based indices, intersected with free
/// flags; order follows the list (duplicates dropped). Out-of-range and
/// fixed indices are skipped.
VectorXi resolveMobileAtoms(const Matter *matter, const std::string &atomList);

/// Explicit candidate indices, intersected with free flags (order preserved,
/// duplicates dropped). Used by pyeonclient compute(..., atoms=...).
VectorXi resolveMobileAtoms(const Matter *matter, const VectorXi &candidates);

/// Pack full (n_atoms,3) rows of mobile atoms into a 3*n_mobile vector.
VectorXd packMobileRows(const AtomMatrix &full, const VectorXi &mobile);

/// Write a 3*n_mobile vector into full AtomMatrix rows (other rows unchanged).
void unpackMobileRows(const VectorXd &packed, const VectorXi &mobile,
                      AtomMatrix &full);

/// Force components on mobile atoms after Matter has a valid force cache
/// (calls getForces under the hood if needed via the non-const Matter API).
VectorXd mobileForces(Matter *matter, const VectorXi &mobile);

} // namespace eonc
