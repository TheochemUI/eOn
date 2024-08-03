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

#include "client/io/ConWriter.hpp"
#include "client/Element.hpp"

namespace eonc::io {

bool ConWriter::writeImpl(const Matter &mat, std::ofstream &fout) {
  // TODO(rg) :: These should be lambda applied sometime
  // if (mat.usePeriodicBoundaries) {
  //     mat.applyPeriodicBoundary();
  // }

  const std::vector<std::string> preHeader{"Generated by eOn\n", "\n"};
  const std::vector<std::string> postHeader{"\n", "\n"};
  size_t numAtoms = mat.numberOfAtoms();
  if (numAtoms == 0) {
    return false;
  }

  std::vector<size_t> componentStartIndices;
  std::vector<double> masses;
  std::vector<size_t> atomicNumbers;

  componentStartIndices.reserve(numAtoms /
                                2); // Estimate to reduce reallocations
  masses.reserve(numAtoms / 2);
  atomicNumbers.reserve(numAtoms / 2);

  componentStartIndices.push_back(0);
  masses.push_back(mat.getMass(0));
  atomicNumbers.push_back(mat.getAtomicNr(0));

  for (size_t i = 1, currentComponent = 0; i < numAtoms; ++i) {
    if (mat.getAtomicNr(i) != atomicNumbers[currentComponent]) {
      ++currentComponent;
      masses.push_back(mat.getMass(i));
      atomicNumbers.push_back(mat.getAtomicNr(i));
      componentStartIndices.push_back(i);
    }
  }
  componentStartIndices.push_back(numAtoms);
  size_t numComponents = masses.size();

  fout << preHeader[0] << preHeader[1];

  Eigen::Vector3d lengths = mat.cell.rowwise().norm();
  fout << fmt::format("{}\t{}\t{}\n", lengths(0), lengths(1), lengths(2));

  Eigen::Vector3d angles;
  angles(0) = std::acos(mat.cell.row(0).dot(mat.cell.row(1)) /
                        (lengths(0) * lengths(1))) *
              180 / M_PI;
  angles(1) = std::acos(mat.cell.row(0).dot(mat.cell.row(2)) /
                        (lengths(0) * lengths(2))) *
              180 / M_PI;
  angles(2) = std::acos(mat.cell.row(1).dot(mat.cell.row(2)) /
                        (lengths(1) * lengths(2))) *
              180 / M_PI;
  fout << fmt::format("{}\t{}\t{}\n", angles(0), angles(1), angles(2));

  fout << postHeader[0] << postHeader[1];

  fout << numComponents << "\n";
  for (size_t j = 0; j < numComponents; ++j) {
    fout << (componentStartIndices[j + 1] - componentStartIndices[j]) << " ";
  }
  fout << "\n";
  for (double mass : masses) {
    fout << mass << " ";
  }
  fout << "\n";
  for (size_t j = 0; j < numComponents; ++j) {
    fout << atomicNumber2symbol(atomicNumbers[j]) << "\n";
    fout << fmt::format("Coordinates of Component {}\n", j + 1);
    for (size_t i = componentStartIndices[j]; i < componentStartIndices[j + 1];
         ++i) {
      fout << fmt::format("{:.17f} {:.17f} {:.17f} {} {:4}\n",
                          mat.getPosition(i, 0), mat.getPosition(i, 1),
                          mat.getPosition(i, 2), mat.getFixed(i), i);
    }
  }

  return true;
}

} // namespace eonc::io
