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

namespace eonc::mat {

void conCell(Matter &mat, const std::array<double, 3> &_lengths,
             const std::array<double, 3> &_degrees) {
  Matrix3S cell = Matrix3S::Zero();
  bool isOrthogonal = std::all_of(_degrees.begin(), _degrees.end(),
                                  [](double angle) { return angle == 90.0; });

  if (isOrthogonal) {
    for (int i = 0; i < 3; ++i) {
      cell(i, i) = _lengths[i];
    }
  } else {
    std::array<double, 3> radians;
    std::transform(_degrees.begin(), _degrees.end(), radians.begin(),
                   [](double angle) { return angle * M_PI / 180.0; });
    // Perform trigonometric calculations on angles to fill cell matrix
    cell(0, 0) = 1.0;
    cell(1, 0) = cos(radians[0]);
    cell(1, 1) = sin(radians[0]);
    cell(2, 0) = cos(radians[1]);
    cell(2, 1) = (cos(radians[2]) - cell(1, 0) * cell(2, 0)) / cell(1, 1);
    cell(2, 2) = sqrt(1.0 - pow(cell(2, 0), 2) - pow(cell(2, 1), 2));

    cell(0, 0) *= _lengths[0];
    cell(1, 0) *= _lengths[1];
    cell(1, 1) *= _lengths[1];
    cell(2, 0) *= _lengths[2];
    cell(2, 1) *= _lengths[2];
    cell(2, 2) *= _lengths[2];
  }

  mat.cell = cell;
  mat.cellInverse = cell.inverse();
  return;
}

void from_con(Matter &mat, const std::string &fname) {
  std::ifstream file(fname);
  if (!file.is_open()) {
    throw std::runtime_error("Failed to open file: " + fname);
  }

  auto readLine = [&file]() -> std::optional<std::string> {
    std::string line;
    if (std::getline(file, line)) {
      return line;
    }
    return std::nullopt;
  };

  auto readArray = [](const std::string &line) -> std::array<double, 3> {
    std::array<double, 3> arr;
    std::istringstream iss(line);
    iss >> arr[0] >> arr[1] >> arr[2];
    return arr;
  };

  std::string headerCon1 = *readLine();
  std::string headerCon2 = *readLine();
  auto lengths = readArray(*readLine());
  auto angles = readArray(*readLine());
  std::string headerCon4 = *readLine();
  std::string headerCon5 = *readLine();

  conCell(mat, lengths, angles);

  size_t Ncomponent{0};
  if (auto line = readLine()) {
    std::istringstream iss(*line);
    if (!(iss >> Ncomponent)) {
      SPDLOG_INFO("The number of components cannot be read. One component is "
                  "assumed instead");
      Ncomponent = 1;
    }
  }

  SPDLOG_INFO("Got {} components", Ncomponent);

  if (Ncomponent < 1) {
    throw std::runtime_error("Unsupported number of components");
  }

  std::vector<size_t> first(Ncomponent + 1, 0);
  if (auto line = readLine()) {
    std::istringstream iss(*line);
    for (size_t j = 0; j < Ncomponent; ++j) {
      size_t Natoms;
      if (!(iss >> Natoms)) {
        throw std::runtime_error("Invalid component count");
      }
      first[j + 1] = Natoms + first[j];
    }
  }

  mat.resize(first[Ncomponent]);

  std::vector<double> mass(Ncomponent);
  if (auto line = readLine()) {
    std::istringstream iss(*line);
    for (size_t j = 0; j < Ncomponent; ++j) {
      if (!(iss >> mass[j])) {
        throw std::runtime_error("Invalid mass entry");
      }
    }
  }

  for (size_t j = 0; j < Ncomponent; ++j) {
    if (auto line = readLine()) {
      std::istringstream iss(*line);
      std::string symbol;
      iss >> symbol;
      int atomicNr = symbol2atomicNumber(symbol);

      // Skip:
      // Coordinates of component blah
      readLine();

      for (size_t i = first[j]; i < first[j + 1]; ++i) {
        mat.setMass(i, mass[j]);
        mat.setAtomicNr(i, atomicNr);
        if (auto line = readLine()) {
          std::istringstream iss(*line);
          if (!(iss >> mat.positions(i, 0) >> mat.positions(i, 1) >>
                mat.positions(i, 2) >> mat.isFixed[i])) {
            throw std::runtime_error("Error parsing atom data");
          }
        }
      }
    }
  }

  if (mat.usePeriodicBoundaries) {
    mat.applyPeriodicBoundary();
  }

  mat.recomputePotential = true;
}

} // namespace eonc::mat
