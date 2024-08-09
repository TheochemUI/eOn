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
#include <fstream>

namespace eonc::mat {

std::optional<std::string> readLine(std::ifstream &file) {
  std::string line;
  if (std::getline(file, line)) {
    return line;
  }
  return std::nullopt;
}

std::array<double, 3> readArray(const std::string &line) {
  std::array<double, 3> arr;
  std::istringstream iss(line);
  iss >> arr[0] >> arr[1] >> arr[2];
  return arr;
}

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

void ConFileParser::con_headers(std::ifstream &file, Matter &mat) {
  headers[0] = *readLine(file);
  headers[1] = *readLine(file);
  auto lengths = readArray(*readLine(file));
  auto angles = readArray(*readLine(file));
  headers[2] = *readLine(file);
  headers[3] = *readLine(file);

  conCell(mat, lengths, angles);

  if (auto line = readLine(file)) {
    std::istringstream iss(*line);
    if (!(iss >> Ncomponent)) {
      SPDLOG_INFO("The number of components cannot be read. One component is "
                  "assumed instead");
      Ncomponent = 1;
    }
  }

  if (Ncomponent < 1) {
    throw std::runtime_error("Unsupported number of components");
  }

  first.resize(Ncomponent + 1, 0);

  size_t totalAtoms{0};
  if (auto line = readLine(file)) {
    std::istringstream iss(*line);
    for (size_t j = 0; j < Ncomponent; ++j) {
      size_t Natoms;
      if (!(iss >> Natoms)) {
        throw std::runtime_error("Invalid component count");
      }
      first[j + 1] = Natoms + first[j];
      totalAtoms += Natoms;
    }
  }

  mat.resize(totalAtoms);
  std::vector<double> mass(totalAtoms);

  if (auto line = readLine(file)) {
    std::istringstream iss(*line);
    for (size_t j = 0; j < Ncomponent; ++j) {
      double componentMass;
      if (!(iss >> componentMass)) {
        throw std::runtime_error("Invalid mass entry");
      }
      std::fill(mass.begin() + first[j], mass.begin() + first[j + 1],
                componentMass);
    }
  }
  mat.setMasses(VectorType::Map(mass.data(), totalAtoms));
}

void ConFileParser::read_cartesian_fix_index(std::ifstream &file, Matter &mat,
                                             AtomMatrix &mdat,
                                             const bool set_atmnr_) {
  int atomicNr = -1;
  for (size_t j = 0; j < Ncomponent; ++j) {
    if (auto line = readLine(file)) {
      if (set_atmnr_) {
        std::istringstream iss(*line);
        std::string symbol;
        iss >> symbol;
        atomicNr = symbol2atomicNumber(symbol);
      }

      // Skip:
      // Coordinates of component blah
      readLine(file);

      for (size_t i = first[j]; i < first[j + 1]; ++i) {
        if (set_atmnr_) {
          mat.setAtomicNr(i, atomicNr);
        }
        if (auto line = readLine(file)) {
          std::istringstream iss(*line);
          if (!(iss >> mdat(i, 0) >> mdat(i, 1) >> mdat(i, 2) >>
                mat.isFixed[i] >> mat.atmindices[i])) {
            throw std::runtime_error("Error parsing data: " + *line);
          }
        }
      }
    }
  }
}

void ConFileParser::finalize_matter(Matter &mat) {
  if (mat.usePeriodicBoundaries) {
    mat.applyPeriodicBoundary();
  }
}

void ConFileParser::parse(Matter &mat, const std::string &fname) {
  std::ifstream _file(fname);
  if (!_file.is_open()) {
    throw std::runtime_error("Failed to open file: " + fname);
  }
  int pos = fname.find_last_of('.');
  if (fname.compare(pos + 1, 3, "con")) {
    con_headers(_file, mat);
    read_cartesian_fix_index(_file, mat, mat.positions, true);
    finalize_matter(mat);
  } else if (fname.compare(pos + 1, 3, "convel")) {
    con_headers(_file, mat);
    read_cartesian_fix_index(_file, mat, mat.positions, true);
    // Empty line between velocities and positions
    readLine(_file);
    read_cartesian_fix_index(_file, mat, mat.velocities, false);
    finalize_matter(mat);
  } else {
    throw std::runtime_error("File extension is not convel or con.");
  }
}

} // namespace eonc::mat
