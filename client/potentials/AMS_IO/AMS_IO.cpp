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

#include "eon/potentials/AMS_IO/AMS_IO.h"
#include "eon/potentials/ExternalCommand.h"

#include <cstddef>
#include <cstring>
#include <format>
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>
#ifndef _WIN32
#include <sys/stat.h>
#include <unistd.h>
#endif

namespace {

// The driver script and its output, both relative to the working directory.
constexpr const char *kRunScript = "run_AMS_IO.sh";
constexpr const char *kOutputFile = "ams_output";
// Truncating rather than appending keeps the previous call's block out of the
// parse, so a failed run cannot hand back stale forces.
constexpr const char *kRunCommand = "./run_AMS_IO.sh > ams_output";

// Lines that precede the two blocks of interest in the AMS output.
constexpr const char *kEnergyMarker = "     CALCULATION RESULTS";
constexpr const char *kGradientMarker =
    "  Index   Atom            d/dx            d/dy            d/dz";

constexpr double kHartreeToEv = 27.2114;
constexpr double kHartreeBohrToEvAngstrom = 51.4220862;

} // namespace

AMS_IO::AMS_IO(const Parameters &p)
    : Potential(PotType::AMS_IO, p) {
  engine = p.ams_options.engine;
  forcefield = p.ams_options.forcefield;
  model = p.ams_options.model;
  xc = p.ams_options.xc;
  return;
}

void AMS_IO::cleanMemory(void) { return; }

AMS_IO::~AMS_IO() { cleanMemory(); }

namespace {

const char *elementArray[] = {
    "Unknown", "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na",
    "Mg",      "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",
    "Cr",      "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br",
    "Kr",      "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag",
    "Cd",      "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr",
    "Nd",      "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
    "Hf",      "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi",
    "Po",      "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",  NULL};

// guess the atom type from the atomic mass,
std::string mass2atom(double atomicmass) {
  return elementArray[int(atomicmass + .5)];
}

int symbol2atomicNumber(char const *symbol) {
  int i = 0;

  while (elementArray[i] != NULL) {
    if (strcmp(symbol, elementArray[i]) == 0) {
      return i;
    }
    i++;
  }
  // invalid symbol
  return -1;
}

char const *atomicNumber2symbol(int n) {
  // The trailing NULL terminates the table, so it bounds the valid range.
  if (n < 0 || static_cast<std::size_t>(n) + 1 >= std::size(elementArray)) {
    throw std::runtime_error(
        std::format("AMS_IO knows no element symbol for atomic number {}", n));
  }
  return elementArray[n];
}
} // namespace

void AMS_IO::force(long N, const double *R, const int *atomicNrs, double *F,
                   double *U, double *variance, const double *box) {
  variance = nullptr;
  passToSystem(N, R, atomicNrs, box);
  // Run a single point AMS_IO calculation and write the results into
  // ams_output
  eonc::pot::runOrThrow(kRunCommand);
  recieveFromSystem(N, F, U);
  return;
}

void AMS_IO::passToSystem(long N, const double *R, const int *atomicNrs,
                          const double *box)
// Creating the standard input file that will be read by the AMS_IO driver
{
  std::ofstream out(kRunScript, std::ios::trunc);
  if (!out) {
    throw std::runtime_error(
        std::format("Could not open {} for writing", kRunScript));
  }

  out << "#!/bin/sh\n";
  out << "ams --delete-old-results <<eor\n";
  out << "Task SinglePoint\n";
  out << "System\n";
  out << " Atoms\n";
  for (long i = 0; i < N; i++) {
    out << std::format("  {}\t{:.19f}\t{:.19f}\t{:.19f}\n",
                       atomicNumber2symbol(atomicNrs[i]), R[i * 3 + 0],
                       R[i * 3 + 1], R[i * 3 + 2]);
  }
  out << " End\n";
  if (!model.empty() || !forcefield.empty()) {
    out << " Lattice\n";
    for (int i = 0; i < 3; i++) {
      out << std::format("  {:.19f}\t{:.19f}\t{:.19f}\n", box[i * 3 + 0],
                         box[i * 3 + 1], box[i * 3 + 2]);
    }
    out << " End\n";
  }
  out << "End\n";
  out << std::format("Engine {}\n", engine);
  if (!forcefield.empty()) {
    out << std::format("     Forcefield {}\n", forcefield);
  }
  if (!model.empty()) {
    out << std::format("     Model {}\n", model);
  }
  if (!xc.empty()) {
    out << std::format("xc {}\n", xc);
    // basis set not specified (default = DZ)
    out << std::format("     hybrid {}\n", xc);
    out << "end\n";
  }
  out << "EndEngine\n";
  out << "Properties\n";
  out << " Gradients\n";
  out << "End\n";
  out << "eor";

  out.close();
  if (!out) {
    throw std::runtime_error(
        std::format("Could not write the AMS input to {}", kRunScript));
  }

#ifndef _WIN32
  if (chmod(kRunScript, S_IRWXU) != 0) {
    throw std::runtime_error(
        std::format("Could not make {} executable", kRunScript));
  }
#endif
  return;
}

void AMS_IO::recieveFromSystem(long N, double *F, double *U) {
  std::ifstream in(kOutputFile);
  if (!in) {
    throw std::runtime_error(
        std::format("Could not open {}; AMS left no output", kOutputFile));
  }

  bool haveEnergy = false;
  bool haveGradients = false;
  std::string line;

  while (std::getline(in, line)) {

    if (line == kEnergyMarker) { // Finding the Energy in the output file
      std::string junk;
      if (!(in >> junk >> junk >> junk >> *U)) {
        throw std::runtime_error(
            std::format("Could not read the energy following \"{}\" in {}",
                        kEnergyMarker, kOutputFile));
      }
      *U = *U * kHartreeToEv; // Energy in hartree to eV
      haveEnergy = true;
    }

    if (line == kGradientMarker) { // Finding the forces
      double index;
      std::string symbol;
      for (long i = 0; i < N; i++) {
        if (!(in >> index >> symbol >> F[i * 3 + 0] >> F[i * 3 + 1] >>
              F[i * 3 + 2])) {
          throw std::runtime_error(
              std::format("{} holds gradients for {} atoms, expected {}",
                          kOutputFile, i, N));
        }
        // AMS_IO gives gradients, not forces, hence the change.
        F[i * 3 + 0] = -F[i * 3 + 0];
        F[i * 3 + 1] = -F[i * 3 + 1];
        F[i * 3 + 2] = -F[i * 3 + 2];
      }
      haveGradients = true;
    }
  }

  if (!haveEnergy || !haveGradients) {
    throw std::runtime_error(std::format(
        "{} holds no {}", kOutputFile,
        !haveEnergy
            ? (!haveGradients ? "energy and no gradient block" : "energy")
            : "gradient block"));
  }

  for (long i = 0; i < 3 * N; i++) {
    F[i] = F[i] * kHartreeBohrToEvAngstrom; // Forces from hartree/bohr to
                                            // eV/Angstrom
  }
  return;
}
