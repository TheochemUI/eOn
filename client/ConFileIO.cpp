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

char const *atomicNumber2symbol(int n) { return elementArray[n]; }

// Strip trailing newlines/carriage returns for readcon header round-tripping
std::string strip_nl(const std::string &s) {
  std::string str(s);
  while (!str.empty() && (str.back() == '\n' || str.back() == '\r'))
    str.pop_back();
  return str;
}

} // namespace

namespace eonc::io {

std::pair<std::array<double, 3>, std::array<double, 3>>
cell_to_lengths_angles(const Matter &m) {
  std::array<double, 3> lengths;
  lengths[0] = m.cell.row(0).norm();
  lengths[1] = m.cell.row(1).norm();
  lengths[2] = m.cell.row(2).norm();
  std::array<double, 3> angles;
  angles[0] =
      eonc::safemath::safe_acos(eonc::safemath::safe_div(
          m.cell.row(0).dot(m.cell.row(1)), lengths[0] * lengths[1])) *
      180.0 / eonc::helpers::pi;
  angles[1] =
      eonc::safemath::safe_acos(eonc::safemath::safe_div(
          m.cell.row(0).dot(m.cell.row(2)), lengths[0] * lengths[2])) *
      180.0 / eonc::helpers::pi;
  angles[2] =
      eonc::safemath::safe_acos(eonc::safemath::safe_div(
          m.cell.row(1).dot(m.cell.row(2)), lengths[1] * lengths[2])) *
      180.0 / eonc::helpers::pi;
  return {lengths, angles};
}

bool matter2con(Matter &m, std::string filename, bool append) {
  int pos = filename.find_last_of('.');
  if (filename.compare(pos + 1, 3, "con")) {
    filename += ".con";
  }

  if (m.usePeriodicBoundaries) {
    m.applyPeriodicBoundary();
  }

  auto [lengths, angles_deg] = cell_to_lengths_angles(m);

  readcon::ConFrameBuilder builder(
      {lengths[0], lengths[1], lengths[2]},
      {angles_deg[0], angles_deg[1], angles_deg[2]},
      {strip_nl(m.headerCon[0]), strip_nl(m.headerCon[1])},
      {strip_nl(m.headerCon[3]), strip_nl(m.headerCon[4])});

  for (long i = 0; i < m.numberOfAtoms(); i++) {
    builder.add_atom(atomicNumber2symbol(m.getAtomicNr(i)), m.getPosition(i, 0),
                     m.getPosition(i, 1), m.getPosition(i, 2),
                     m.getFixed(i) != 0, static_cast<uint64_t>(m.atomIndex(i)),
                     m.getMass(i));
  }

  auto frame = builder.build();
  std::vector<readcon::ConFrame> frames;
  if (append) {
    try {
      frames = readcon::read_all_frames(filename);
    } catch (...) {
      // File doesn't exist yet, start fresh
    }
  }
  frames.push_back(std::move(frame));
  readcon::ConFrameWriter writer(filename, 17);
  writer.extend(frames);
  return true;
}

// Load atomic coordinates from a .con file via readcon-core (mmap reader)
bool con2matter(Matter &m, std::string filename) {
  int pos = filename.find_last_of('.');
  if (filename.compare(pos + 1, 3, "con")) {
    filename += ".con";
  }
  try {
    auto frame = readcon::read_first_frame(filename);
    return con2matter(m, frame);
  } catch (const std::exception &e) {
    EONC_LOG_ERROR("Failed to read {}: {}", filename, e.what());
    return false;
  }
}

// Populate Matter from a parsed readcon frame
bool con2matter(Matter &m, const readcon::ConFrame &frame) {
  const auto &atoms = frame.atoms();
  const auto &lengths = frame.cell();
  const auto &angles_deg = frame.angles();
  const auto &prebox = frame.prebox_header();
  const auto &postbox = frame.postbox_header();

  // Store headers for round-tripping via matter2con
  m.headerCon[0] = prebox[0] + "\n";
  m.headerCon[1] = prebox[1] + "\n";
  m.headerCon[3] = postbox[0] + "\n";
  m.headerCon[4] = postbox[1] + "\n";

  // Build cell matrix from lengths and angles
  double angles[3] = {angles_deg[0], angles_deg[1], angles_deg[2]};
  if (angles[0] == 90.0 && angles[1] == 90.0 && angles[2] == 90.0) {
    m.cell.setZero();
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

  // Store angles line for convel round-tripping
  m.headerCon[2] =
      std::format("{} {} {}\n", angles_deg[0], angles_deg[1], angles_deg[2]);

  m.resize(static_cast<long>(atoms.size()));
  for (size_t i = 0; i < atoms.size(); ++i) {
    m.positions(i, 0) = atoms[i].x;
    m.positions(i, 1) = atoms[i].y;
    m.positions(i, 2) = atoms[i].z;
    m.setMass(i, atoms[i].mass);
    m.setAtomicNr(i, static_cast<int>(atoms[i].atomic_number));
    m.setFixed(i, atoms[i].is_fixed ? 1 : 0);
    m.atomIndex(i) = static_cast<int>(atoms[i].atom_id);
  }

  if (m.usePeriodicBoundaries) {
    m.applyPeriodicBoundary();
  }
  m.recomputePotential = true;
  return true;
}

bool matter2convel(Matter &m, std::string filename) {
  int pos = filename.find_last_of('.');
  if (filename.compare(pos + 1, 6, "convel")) {
    filename += ".convel";
  }

  if (m.usePeriodicBoundaries) {
    m.applyPeriodicBoundary();
  }

  auto [lengths, angles_deg] = cell_to_lengths_angles(m);

  readcon::ConFrameBuilder builder(
      {lengths[0], lengths[1], lengths[2]},
      {angles_deg[0], angles_deg[1], angles_deg[2]},
      {strip_nl(m.headerCon[0]), strip_nl(m.headerCon[1])},
      {strip_nl(m.headerCon[3]), strip_nl(m.headerCon[4])});

  for (long i = 0; i < m.numberOfAtoms(); i++) {
    builder.add_atom_with_velocity(
        atomicNumber2symbol(m.getAtomicNr(i)), m.getPosition(i, 0),
        m.getPosition(i, 1), m.getPosition(i, 2), m.getFixed(i) != 0,
        static_cast<uint64_t>(m.atomIndex(i)), m.getMass(i),
        m.velocities(i, 0), m.velocities(i, 1), m.velocities(i, 2));
  }

  auto frame = builder.build();
  readcon::ConFrameWriter writer(filename, 6);
  std::vector<readcon::ConFrame> frames;
  frames.push_back(std::move(frame));
  writer.extend(frames);
  return true;
}

bool convel2matter(Matter &m, std::string filename) {
  int pos = filename.find_last_of('.');
  if (filename.compare(pos + 1, 6, "convel")) {
    filename += ".convel";
  }
  try {
    auto frame = readcon::read_first_frame(filename);
    bool ok = con2matter(m, frame);
    if (!ok)
      return false;
    if (frame.has_velocities()) {
      const auto &atoms = frame.atoms();
      for (size_t i = 0; i < atoms.size(); ++i) {
        m.setVelocity(i, 0, atoms[i].vx);
        m.setVelocity(i, 1, atoms[i].vy);
        m.setVelocity(i, 2, atoms[i].vz);
      }
    }
    return true;
  } catch (const std::exception &e) {
    EONC_LOG_ERROR("Failed to read convel {}: {}", filename, e.what());
    return false;
  }
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
