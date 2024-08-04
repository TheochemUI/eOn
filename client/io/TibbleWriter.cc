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

#include "client/io/TibbleWriter.hpp"
#include "client/Element.hpp"

namespace eonc::io {

/**
 * @brief Writes the Matter object to a file suitable for reading in as a Tibble
 * in R.
 *
 * This function writes the details of the Matter object, such as atomic
 * positions, forces, and atomic types, to a file.
 *
 * @param mat The Matter object containing atomic data.
 * @param fout The output file stream where the data will be written.
 * @return True if the file is written successfully, false otherwise.
 *
 * @note The function assumes that atoms are grouped by their atomic numbers.
 *       If the Matter object is empty (i.e., contains no atoms), the function
 *       returns false.
 */
bool TibbleWriter::writeImpl(const Matter &mat, std::ofstream &fout) {
  // TODO(rg) :: These should be lambda applied sometime
  // if (mat.usePeriodicBoundaries) {
  //     mat.applyPeriodicBoundary();
  // }
  using namespace fmt::literals;
  const size_t nAtoms = mat.numberOfAtoms();
  const AtomMatrix fSys = const_cast<Matter &>(mat).getForces();
  const double eSys = const_cast<Matter &>(mat).getPotentialEnergy();
  const AtomMatrix pos = mat.getPositions();
  fout << "x y z fx fy fz energy mass symbol atmID fixed\n";
  for (size_t idx{0}; idx < nAtoms; idx++) {
    fout << fmt::format(
        "{x} {y} {z} {fx} {fy} {fz} {energy} {mass} {symbol} {idx} {fixed}\n",
        "x"_a = pos.row(idx)[0], "y"_a = pos.row(idx)[1],
        "z"_a = pos.row(idx)[2], "fx"_a = fSys.row(idx)[0],
        "fy"_a = fSys.row(idx)[1], "fz"_a = fSys.row(idx)[2], "energy"_a = eSys,
        "mass"_a = mat.getMass(idx),
        "symbol"_a = atomicNumber2symbol(mat.getAtomicNr(idx)),
        "idx"_a = (idx + 1),
        /* NOTE(rg): idx MAY not be the same id as before */
        "fixed"_a = mat.getFixed(idx));
  }

  return true;
}

} // namespace eonc::io
