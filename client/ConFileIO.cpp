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
#include "ConFileIO.h"
#include "EonLogger.h"
#include "HelperFunctions.h"
#include "Matter.h"
#include "SafeMath.h"

#include <cstring>
#include <format>
#include <fstream>
#include <iostream>
#include <string>

namespace {

const char *elementArray[] = {
    "Unknown", "H",  "He", "Li", "Be", "B",    "C",  "N",  "O",  "F",  "Ne",
    "Na",      "Mg", "Al", "Si", "P",  "S",    "Cl", "Ar", "K",  "Ca", "Sc",
    "Ti",      "V",  "Cr", "Mn", "Fe", "Co",   "Ni", "Cu", "Zn", "Ga", "Ge",
    "As",      "Se", "Br", "Kr", "Rb", "Sr",   "Y",  "Zr", "Nb", "Mo", "Tc",
    "Ru",      "Rh", "Pd", "Ag", "Cd", "In",   "Sn", "Sb", "Te", "I",  "Xe",
    "Cs",      "Ba", "La", "Ce", "Pr", "Nd",   "Pm", "Sm", "Eu", "Gd", "Tb",
    "Dy",      "Ho", "Er", "Tm", "Yb", "Lu",   "Hf", "Ta", "W",  "Re", "Os",
    "Ir",      "Pt", "Au", "Hg", "Tl", "Pb",   "Bi", "Po", "At", "Rn", "Fr",
    "Ra",      "Ac", "Th", "Pa", "U",  nullptr};

const int MAXC = 100;

int symbol2atomicNumber(char const *symbol) {
  int i = 0;
  while (elementArray[i] != nullptr) {
    if (strcmp(symbol, elementArray[i]) == 0) {
      return i;
    }
    i++;
  }
  return -1;
}

char const *atomicNumber2symbol(int n) { return elementArray[n]; }

} // namespace

namespace eonc::io {

bool matter2con(Matter &m, std::string filename, bool append) {
  FILE *file;
  int pos = filename.find_last_of('.');
  if (filename.compare(pos + 1, 3, "con")) {
    filename += ".con";
  };
  if (append) {
    file = fopen(filename.c_str(), "ab");
  } else {
    file = fopen(filename.c_str(), "wb");
  }
  bool state = matter2con(m, file);
  fclose(file);
  return state;
}

bool matter2con(Matter &m, FILE *file) {
  long int i;
  int j;
  long int Nfix = 0;
  int Ncomponent = 0;
  int first[MAXC];
  double mass[MAXC];
  long atomicNrs[MAXC];
  first[0] = 0;

  if (m.usePeriodicBoundaries) {
    m.applyPeriodicBoundary();
  }

  if (m.numberOfAtoms() > 0) {
    if (m.getFixed(0))
      Nfix = 1;
    mass[0] = m.getMass(0);
    atomicNrs[0] = m.getAtomicNr(0);
  };
  j = 0;
  for (i = 1; i < m.numberOfAtoms(); i++) {
    if (m.getFixed(i))
      Nfix++;
    if (m.getAtomicNr(i) != atomicNrs[j]) {
      j++;
      if (j >= MAXC) {
        EONC_LOG_ERROR("Does not support more than {} components and the "
                       "atoms must be ordered by component.",
                       MAXC);
        return false;
      };
      mass[j] = m.getMass(i);
      atomicNrs[j] = m.getAtomicNr(i);
      first[j] = i;
    }
  }
  first[j + 1] = m.numberOfAtoms();
  Ncomponent = j + 1;

  fputs(m.headerCon[0].c_str(), file);
  fputs(m.headerCon[1].c_str(), file);
  double lengths[3];
  lengths[0] = m.cell.row(0).norm();
  lengths[1] = m.cell.row(1).norm();
  lengths[2] = m.cell.row(2).norm();
  fprintf(file, "%f\t%f\t%f\n", lengths[0], lengths[1], lengths[2]);
  double angles[3];
  angles[0] = eonc::safemath::safe_acos(eonc::safemath::safe_div(
                  m.cell.row(0).dot(m.cell.row(1)), lengths[0] * lengths[1])) *
              180 / eonc::helpers::pi;
  angles[1] = eonc::safemath::safe_acos(eonc::safemath::safe_div(
                  m.cell.row(0).dot(m.cell.row(2)), lengths[0] * lengths[2])) *
              180 / eonc::helpers::pi;
  angles[2] = eonc::safemath::safe_acos(eonc::safemath::safe_div(
                  m.cell.row(1).dot(m.cell.row(2)), lengths[1] * lengths[2])) *
              180 / eonc::helpers::pi;
  fprintf(file, "%f\t%f\t%f\n", angles[0], angles[1], angles[2]);
  fputs(m.headerCon[3].c_str(), file);
  fputs(m.headerCon[4].c_str(), file);

  fprintf(file, "%d\n", Ncomponent);
  for (j = 0; j < Ncomponent; j++) {
    fprintf(file, "%d ", first[j + 1] - first[j]);
  }
  fprintf(file, "\n");
  for (j = 0; j < Ncomponent; j++) {
    fprintf(file, "%f ", mass[j]);
  }
  fprintf(file, "\n");
  for (j = 0; j < Ncomponent; j++) {
    fprintf(file, "%s\n", atomicNumber2symbol(atomicNrs[j]));
    fprintf(file, "Coordinates of Component %d\n", j + 1);
    for (i = first[j]; i < first[j + 1]; i++) {
      fprintf(file, "%22.17f %22.17f %22.17f %d %4d\n", m.getPosition(i, 0),
              m.getPosition(i, 1), m.getPosition(i, 2), m.getFixed(i),
              m.atomIndex(i));
    }
  }
  return true;
}

bool con2matter(Matter &m, std::string filename) {
  FILE *file;
  int pos = filename.find_last_of('.');
  if (filename.compare(pos + 1, 3, "con")) {
    filename += ".con";
  }
  file = fopen(filename.c_str(), "rb");
  if (!file) {
    EONC_LOG_ERROR("File {} was not found.", filename);
    return false;
  }
  bool state = con2matter(m, file);
  fclose(file);
  return state;
}

bool con2matter(Matter &m, FILE *file) {
  char line[255];
  fgets(line, sizeof(line), file);
  m.headerCon[0] = line;

  long int i;
  int j;

  fgets(line, sizeof(line), file);
  m.headerCon[1] = line;

  double lengths[3];
  fgets(line, sizeof(line), file);
  sscanf(line, "%lf %lf %lf", &lengths[0], &lengths[1], &lengths[2]);

  double angles[3];
  fgets(line, sizeof(line), file);
  m.headerCon[2] = line;
  sscanf(line, "%lf %lf %lf", &angles[0], &angles[1], &angles[2]);

  if (angles[0] == 90.0 && angles[1] == 90.0 && angles[2] == 90.0) {
    m.cell(0, 0) = lengths[0];
    m.cell(1, 1) = lengths[1];
    m.cell(2, 2) = lengths[2];
  } else {
    angles[0] *= eonc::helpers::pi / 180.0;
    angles[1] *= eonc::helpers::pi / 180.0;
    angles[2] *= eonc::helpers::pi / 180.0;

    m.cell(0, 0) = 1.0;
    m.cell(1, 0) = cos(angles[0]);
    m.cell(1, 1) = sin(angles[0]);
    m.cell(2, 0) = cos(angles[1]);
    m.cell(2, 1) =
        (cos(angles[2]) - m.cell(1, 0) * m.cell(2, 0)) / m.cell(1, 1);
    m.cell(2, 2) = eonc::safemath::safe_sqrt(1.0 - pow(m.cell(2, 0), 2) -
                                             pow(m.cell(2, 1), 2));

    m.cell(0, 0) *= lengths[0];
    m.cell(1, 0) *= lengths[1];
    m.cell(1, 1) *= lengths[1];
    m.cell(2, 0) *= lengths[2];
    m.cell(2, 1) *= lengths[2];
    m.cell(2, 2) *= lengths[2];
  }
  m.cellInverse = m.cell.inverse();

  fgets(line, sizeof(line), file);
  m.headerCon[3] = line;
  fgets(line, sizeof(line), file);
  m.headerCon[4] = line;

  fgets(line, sizeof(line), file);
  int Ncomponent;
  if (sscanf(line, "%d", &Ncomponent) == 0) {
    EONC_LOG_INFO("The number of components cannot be read. One "
                  "component is assumed instead");
    Ncomponent = 1;
  }
  if ((Ncomponent > MAXC) || (Ncomponent < 1)) {
    EONC_LOG_ERROR(
        "con2atoms doesn't support more that {} components or less than 1",
        MAXC);
    return false;
  }

  long int first[MAXC + 1];
  long int Natoms = 0;
  first[0] = 0;

  fgets(line, sizeof(line), file);
  char *split = strtok(line, " \t");
  for (j = 0; j < Ncomponent; j++) {
    if (split == nullptr) {
      EONC_LOG_ERROR(
          "input con file does not list the number of each component");
      return false;
    }
    if (sscanf(split, "%ld", &Natoms) != 1) {
      EONC_LOG_ERROR(
          "input con file does not list the number of each component");
      return false;
    }
    first[j + 1] = Natoms + first[j];
    split = strtok(nullptr, " \t");
  }

  m.resize(first[Ncomponent]);

  double mass[MAXC];
  fgets(line, sizeof(line), file);
  split = strtok(line, " \t");

  for (j = 0; j < Ncomponent; j++) {
    if (split == nullptr) {
      EONC_LOG_ERROR("input con file does not list enough masses");
      return false;
    }
    if (sscanf(split, "%lf", &mass[j]) != 1) {
      EONC_LOG_ERROR("input con file does not list enough masses");
      return false;
    }
    split = strtok(nullptr, " \t");
  }

  int atomicNr;
  int fixed;
  double x, y, z;
  for (j = 0; j < Ncomponent; j++) {
    char symbol[3];
    fgets(line, sizeof(line), file);
    sscanf(line, "%2s", symbol);
    atomicNr = symbol2atomicNumber(symbol);
    fgets(line, sizeof(line), file); // skip one line
    for (i = first[j]; i < first[j + 1]; i++) {
      m.setMass(i, mass[j]);
      m.setAtomicNr(i, atomicNr);
      fgets(line, sizeof(line), file);
      if (strlen(line) < 6) {
        EONC_LOG_ERROR("error parsing position in con file");
        return false;
      }

      long origIdx = i; // default: sequential
      sscanf(line, "%lf %lf %lf %d %ld", &x, &y, &z, &fixed, &origIdx);
      m.positions(i, 0) = x;
      m.positions(i, 1) = y;
      m.positions(i, 2) = z;
      m.setFixed(i, static_cast<bool>(fixed));
      m.atomIndex(i) = static_cast<int>(origIdx);
    }
  }
  if (m.usePeriodicBoundaries) {
    m.applyPeriodicBoundary();
  }
  m.recomputePotential = true;
  return true;
}

bool matter2convel(Matter &m, std::string filename) {
  FILE *file;
  int pos = filename.find_last_of('.');
  if (filename.compare(pos + 1, 6, "convel")) {
    filename += ".convel";
  }
  file = fopen(filename.c_str(), "w");
  bool state = matter2convel(m, file);
  fclose(file);
  return state;
}

bool matter2convel(Matter &m, FILE *file) {
  long int i;
  int j;
  long int Nfix = 0;
  int Ncomponent = 0;
  int first[MAXC];
  double mass[MAXC];
  long atomicNrs[MAXC];
  first[0] = 0;

  if (m.usePeriodicBoundaries) {
    m.applyPeriodicBoundary();
  }

  if (m.numberOfAtoms() > 0) {
    if (m.getFixed(0))
      Nfix = 1;
    mass[0] = m.getMass(0);
    atomicNrs[0] = m.getAtomicNr(0);
  }
  j = 0;
  for (i = 1; i < m.numberOfAtoms(); i++) {
    if (m.getFixed(i))
      Nfix++;
    if (m.getAtomicNr(i) != atomicNrs[j]) {
      j++;
      if (j >= MAXC) {
        EONC_LOG_ERROR("Does not support more than {} components and the "
                       "atoms must be ordered by component.",
                       MAXC);
        return false;
      }
      mass[j] = m.getMass(i);
      atomicNrs[j] = m.getAtomicNr(i);
      first[j] = i;
    }
  }
  first[j + 1] = m.numberOfAtoms();
  Ncomponent = j + 1;

  fputs(m.headerCon[0].c_str(), file);
  fputs(m.headerCon[1].c_str(), file);
  double lengths[3];
  lengths[0] = m.cell.row(0).norm();
  lengths[1] = m.cell.row(1).norm();
  lengths[2] = m.cell.row(2).norm();
  fprintf(file, "%f\t%f\t%f\n", lengths[0], lengths[1], lengths[2]);
  double angles[3];
  angles[0] = eonc::safemath::safe_acos(eonc::safemath::safe_div(
                  m.cell.row(0).dot(m.cell.row(1)), lengths[0] * lengths[1])) *
              180 / eonc::helpers::pi;
  angles[1] = eonc::safemath::safe_acos(eonc::safemath::safe_div(
                  m.cell.row(0).dot(m.cell.row(2)), lengths[0] * lengths[2])) *
              180 / eonc::helpers::pi;
  angles[2] = eonc::safemath::safe_acos(eonc::safemath::safe_div(
                  m.cell.row(1).dot(m.cell.row(2)), lengths[1] * lengths[2])) *
              180 / eonc::helpers::pi;
  fprintf(file, "%f\t%f\t%f\n", angles[0], angles[1], angles[2]);
  fputs(m.headerCon[3].c_str(), file);
  fputs(m.headerCon[4].c_str(), file);

  fprintf(file, "%d\n", Ncomponent);
  for (j = 0; j < Ncomponent; j++) {
    fprintf(file, "%d ", first[j + 1] - first[j]);
  }
  fprintf(file, "\n");
  for (j = 0; j < Ncomponent; j++) {
    fprintf(file, "%f ", mass[j]);
  }
  fprintf(file, "\n");
  for (j = 0; j < Ncomponent; j++) {
    fprintf(file, "%s\n", atomicNumber2symbol(atomicNrs[j]));
    fprintf(file, "Coordinates of Component %d\n", j + 1);
    for (i = first[j]; i < first[j + 1]; i++) {
      fprintf(file, "%11.6f\t%11.6f\t%11.6f\t%d\t%ld\n", m.getPosition(i, 0),
              m.getPosition(i, 1), m.getPosition(i, 2), m.getFixed(i), i);
    }
  }
  fprintf(file, "\n");
  for (j = 0; j < Ncomponent; j++) {
    fprintf(file, "%s\n", atomicNumber2symbol(atomicNrs[j]));
    fprintf(file, "Velocities of Component %d\n", j + 1);
    for (i = first[j]; i < first[j + 1]; i++) {
      fprintf(file, "%11.6f\t%11.6f\t%11.6f\t%d\t%ld\n", m.velocities(i, 0),
              m.velocities(i, 1), m.velocities(i, 2), m.getFixed(i), i);
    }
  }
  return true;
}

bool convel2matter(Matter &m, std::string filename) {
  FILE *file;
  int pos = filename.find_last_of('.');
  if (filename.compare(pos + 1, 6, "convel")) {
    filename += ".convel";
  }
  file = fopen(filename.c_str(), "rb");
  if (!file) {
    EONC_LOG_ERROR("File {} was not found.", filename);
    return false;
  }
  bool state = convel2matter(m, file);
  fclose(file);
  return state;
}

bool convel2matter(Matter &m, FILE *file) {
  char line[255];
  fgets(line, sizeof(line), file);
  m.headerCon[0] = line;

  long int i;
  int j;

  fgets(line, sizeof(line), file);
  m.headerCon[1] = line;

  double lengths[3];
  fgets(line, sizeof(line), file);
  sscanf(line, "%lf %lf %lf", &lengths[0], &lengths[1], &lengths[2]);

  double angles[3];
  fgets(line, sizeof(line), file);
  m.headerCon[2] = line;
  sscanf(line, "%lf %lf %lf", &angles[0], &angles[1], &angles[2]);

  if (angles[0] == 90.0 && angles[1] == 90.0 && angles[2] == 90.0) {
    m.cell(0, 0) = lengths[0];
    m.cell(1, 1) = lengths[1];
    m.cell(2, 2) = lengths[2];
  } else {
    angles[0] *= eonc::helpers::pi / 180.0;
    angles[1] *= eonc::helpers::pi / 180.0;
    angles[2] *= eonc::helpers::pi / 180.0;

    m.cell(0, 0) = 1.0;
    m.cell(1, 0) = cos(angles[0]);
    m.cell(1, 1) = sin(angles[0]);
    m.cell(2, 0) = cos(angles[1]);
    m.cell(2, 1) =
        (cos(angles[2]) - m.cell(1, 0) * m.cell(2, 0)) / m.cell(1, 1);
    m.cell(2, 2) = eonc::safemath::safe_sqrt(1.0 - pow(m.cell(2, 0), 2) -
                                             pow(m.cell(2, 1), 2));

    m.cell(0, 0) *= lengths[0];
    m.cell(1, 0) *= lengths[1];
    m.cell(1, 1) *= lengths[1];
    m.cell(2, 0) *= lengths[2];
    m.cell(2, 1) *= lengths[2];
    m.cell(2, 2) *= lengths[2];
  }
  m.cellInverse = m.cell.inverse();

  fgets(line, sizeof(line), file);
  m.headerCon[3] = line;
  fgets(line, sizeof(line), file);
  m.headerCon[4] = line;

  fgets(line, sizeof(line), file);
  int Ncomponent;
  if (sscanf(line, "%d", &Ncomponent) == 0) {
    EONC_LOG_INFO("The number of components cannot be read. One "
                  "component is assumed instead");
    Ncomponent = 1;
  }
  if ((Ncomponent > MAXC) || (Ncomponent < 1)) {
    EONC_LOG_ERROR(
        "con2atoms doesn't support more that {} components or less than 1",
        MAXC);
    return false;
  }

  long int first[MAXC + 1];
  long int Natoms = 0;
  first[0] = 0;

  for (j = 0; j < Ncomponent; j++) {
    fscanf(file, "%ld", &Natoms);
    first[j + 1] = Natoms + first[j];
  }

  fgets(line, sizeof(line), file);
  m.resize(first[Ncomponent]);
  double mass[MAXC];
  for (j = 0; j < Ncomponent; j++) {
    fscanf(file, "%lf", &mass[j]);
  }

  fgets(line, sizeof(line), file);
  int atomicNr;
  int fixed;
  double x, y, z;
  for (j = 0; j < Ncomponent; j++) {
    char symbol[3];
    fgets(line, sizeof(line), file);
    sscanf(line, "%2s\n", symbol);
    atomicNr = symbol2atomicNumber(symbol);
    fgets(line, sizeof(line), file); // skip one line
    for (i = first[j]; i < first[j + 1]; i++) {
      m.setMass(i, mass[j]);
      m.setAtomicNr(i, atomicNr);
      fgets(line, sizeof(line), file);
      sscanf(line, "%lf %lf %lf %d\n", &x, &y, &z, &fixed);
      m.setPosition(i, 0, x);
      m.setPosition(i, 1, y);
      m.setPosition(i, 2, z);
      m.setFixed(i, static_cast<bool>(fixed));
    }
  }

  fgets(line, sizeof(line), file);
  for (j = 0; j < Ncomponent; j++) {
    fgets(line, sizeof(line), file);
    fgets(line, sizeof(line), file); // skip one line
    for (i = first[j]; i < first[j + 1]; i++) {
      fgets(line, sizeof(line), file);
      sscanf(line, "%lf %lf %lf %d\n", &x, &y, &z, &fixed);
      m.setVelocity(i, 0, x);
      m.setVelocity(i, 1, y);
      m.setVelocity(i, 2, z);
    }
  }

  if (m.usePeriodicBoundaries) {
    m.applyPeriodicBoundary();
  }
  return true;
}

void matter2xyz(Matter &m, std::string filename, bool append) {
  FILE *file;
  long int i;
  filename += ".xyz";
  if (append) {
    file = fopen(filename.c_str(), "ab");
  } else {
    file = fopen(filename.c_str(), "wb");
  }
  if (file == 0) {
    std::cerr << "Can't create file " << filename << std::endl;
    exit(1);
  }
  fprintf(file, "%ld\nGenerated by eOn\n", m.numberOfAtoms());

  if (m.usePeriodicBoundaries) {
    m.applyPeriodicBoundary();
  }

  for (i = 0; i < m.numberOfAtoms(); i++) {
    fprintf(file, "%s\t%11.6f\t%11.6f\t%11.6f\n",
            atomicNumber2symbol(m.getAtomicNr(i)), m.getPosition(i, 0),
            m.getPosition(i, 1), m.getPosition(i, 2));
  }
  fclose(file);
}

void writeTibble(Matter &m, std::string fname) {
  AtomMatrix fSys = m.getForces();
  std::ofstream out(fname);
  double eSys = m.getPotentialEnergy();
  AtomMatrix pos = m.getPositions();
  out << "x y z fx fy fz energy mass symbol atmID fixed\n";
  for (auto idx{0}; idx < m.numberOfAtoms(); idx++) {
    out << std::format("{} {} {} {} {} {} {} {} {} {} {}\n", pos.row(idx)[0],
                       pos.row(idx)[1], pos.row(idx)[2], fSys.row(idx)[0],
                       fSys.row(idx)[1], fSys.row(idx)[2], eSys, m.getMass(idx),
                       atomicNumber2symbol(m.getAtomicNr(idx)), (idx + 1),
                       m.getFixed(idx));
  }
}

} // namespace eonc::io
