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
#include "client/matter/MatterCreator.hpp"
#include "client/Eigen.h"
#include "client/Element.hpp"
#include "magic_enum/magic_enum.hpp"
#include "readCon/include/ReadCon.hpp"

namespace eonc::mat {

void from_con(Matter &mat, const std::string &fname) {
  std::vector<std::string> fconts =
      yodecon::helpers::file::read_con_file(fname);
  auto singleCon =
      yodecon::create_single_con<yodecon::types::ConFrameVec>(fconts);
  const size_t natoms = singleCon.atom_id.size();
  mat.nAtoms = natoms;
  mat.resize(natoms);

  // Set the cell and its inverse
  // TODO(rg) : Refactor out
  Matrix3S cell{Matrix3S::Zero()}, cellInverse{Matrix3S::Zero()};
  if (singleCon.angles[0] == 90.0 && singleCon.angles[1] == 90.0 &&
      singleCon.angles[2] == 90.0) {
    cell(0, 0) = singleCon.boxl[0];
    cell(1, 1) = singleCon.boxl[1];
    cell(2, 2) = singleCon.boxl[2];
  } else {
    singleCon.angles[0] *= M_PI / 180.0;
    singleCon.angles[1] *= M_PI / 180.0;
    singleCon.angles[2] *= M_PI / 180.0;

    cell(0, 0) = 1.0;
    cell(1, 0) = cos(singleCon.angles[0]);
    cell(1, 1) = sin(singleCon.angles[0]);
    cell(2, 0) = cos(singleCon.angles[1]);
    cell(2, 1) =
        (cos(singleCon.angles[2]) - cell(1, 0) * cell(2, 0)) / cell(1, 1);
    cell(2, 2) = sqrt(1.0 - pow(cell(2, 0), 2) - pow(cell(2, 1), 2));

    cell(0, 0) *= singleCon.boxl[0];
    cell(1, 0) *= singleCon.boxl[1];
    cell(1, 1) *= singleCon.boxl[1];
    cell(2, 0) *= singleCon.boxl[2];
    cell(2, 1) *= singleCon.boxl[2];
    cell(2, 2) *= singleCon.boxl[2];
  }
  cellInverse = cell.inverse();
  mat.cell = cell;
  mat.cellInverse = cellInverse;

  // Set other quantities
  AtomMatrix pos(natoms, 3);
  VectorType masses(natoms);
  Vector<size_t> atmnrs(natoms);
  Vector<int> isFixed(natoms);
  for (size_t idx{0}; idx < natoms; idx++) {
    // Needed since std::vector<bool> doesn't have .data()
    isFixed[idx] = static_cast<int>(singleCon.is_fixed[idx]);
    pos(idx, 0) = singleCon.x[idx];
    pos(idx, 1) = singleCon.y[idx];
    pos(idx, 2) = singleCon.z[idx];
    const auto elem =
        magic_enum::enum_cast<Element>(singleCon.symbol[idx]).value();
    // value() will fail anyway if an unsupported element is found
    const auto edat = elementData.find(elem)->second;
    masses[idx] = edat.atomicMass;
    atmnrs[idx] = edat.atomicNumber;
  }
  mat.setMasses(masses);
  mat.setAtomicNrs(atmnrs);
  mat.setFixedMask(isFixed);
  // Automatically applies periodic boundaries and sets recomputePotential
  mat.setPositions(pos);
}

} // namespace eonc::mat
