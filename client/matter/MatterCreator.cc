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

  char line[255];
  std::string headerCon1, headerCon2, headerCon4, headerCon5, headerCon6;

  file.getline(line, sizeof(line));
  headerCon1 = line;
  file.getline(line, sizeof(line));
  headerCon2 = line;

  std::array<double, 3> lengths;
  file.getline(line, sizeof(line));
  sscanf(line, "%lf %lf %lf", &lengths[0], &lengths[1], &lengths[2]);

  std::array<double, 3> angles;
  file.getline(line, sizeof(line));
  headerCon4 = line;
  sscanf(headerCon4.c_str(), "%lf %lf %lf", &angles[0], &angles[1], &angles[2]);

  conCell(mat, lengths, angles);

  file.getline(line, sizeof(line));
  headerCon5 = line;
  file.getline(line, sizeof(line));
  headerCon6 = line;

  file.getline(line, sizeof(line));
  size_t Ncomponent;
  if (sscanf(line, "%zu", &Ncomponent) != 1) {
    SPDLOG_INFO("The number of components cannot be read. One "
                "component is assumed instead");
    Ncomponent = 1;
  }

  size_t first[101];
  size_t Natoms = 0;
  first[0] = 0;

  file.getline(line, sizeof(line));
  char *split = strtok(line, " \t");
  for (size_t j = 0; j < Ncomponent; j++) {
    if (!split || sscanf(split, "%zu", &Natoms) != 1) {
      throw std::runtime_error("Invalid component count");
    }
    first[j + 1] = Natoms + first[j];
    split = strtok(NULL, " \t");
  }

  mat.resize(first[Ncomponent]);

  double mass[100];
  file.getline(line, sizeof(line));
  split = strtok(line, " \t");
  for (size_t j = 0; j < Ncomponent; j++) {
    if (!split || sscanf(split, "%lf", &mass[j]) != 1) {
      throw std::runtime_error("Invalid mass entry");
    }
    split = strtok(NULL, " \t");
  }

  for (size_t j = 0; j < Ncomponent; j++) {
    char symbol[3];
    file.getline(line, sizeof(line));
    sscanf(line, "%2s", symbol);
    int atomicNr = symbol2atomicNumber(symbol);
    file.getline(line, sizeof(line)); // skip one line
    for (size_t i = first[j]; i < first[j + 1]; i++) {
      mat.setMass(i, mass[j]);
      mat.setAtomicNr(i, atomicNr);
      file.getline(line, sizeof(line));
      if (sscanf(line, "%lf %lf %lf %d", &mat.positions(i, 0),
                 &mat.positions(i, 1), &mat.positions(i, 2),
                 &mat.isFixed[i]) != 4) {
        throw std::runtime_error("Error parsing atom data");
      }
    }
  }

  if (mat.usePeriodicBoundaries) {
    mat.applyPeriodicBoundary();
  }

  mat.recomputePotential = true;
}

} // namespace eonc::mat
