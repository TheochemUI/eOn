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
#include "Hessian.h"
#include "EonLogger.h"
#include "HelperFunctions.h"
#include "SafeMath.h"

#include <cmath>
#include <fstream>
#include <sstream>
#include <string>

namespace {

// atom_list entries are *mobile / displaced* atoms for FD (hybrid/PHVA-class
// active set). Intersect with non-fixed atoms in HessianJob.

bool isCentralScheme(const std::string &scheme) {
  return scheme == "central" || scheme == "CENTRAL" || scheme == "Central";
}

// Checkpoint: first line "eon_hess_ckpt <size> <next_col>", then size*size
// doubles in row-major order matching MatrixXd storage.
bool loadColumnCheckpoint(const std::string &path, int size, int &nextCol,
                          MatrixXd &H) {
  std::ifstream in(path);
  if (!in) {
    return false;
  }
  std::string tag;
  int fileSize = 0;
  in >> tag >> fileSize >> nextCol;
  if (!in || tag != "eon_hess_ckpt" || fileSize != size || nextCol < 0 ||
      nextCol > size) {
    return false;
  }
  H.resize(size, size);
  for (int i = 0; i < size; ++i) {
    for (int j = 0; j < size; ++j) {
      double v = 0.0;
      in >> v;
      if (!in) {
        return false;
      }
      H(i, j) = v;
    }
  }
  return true;
}

bool saveColumnCheckpoint(const std::string &path, int size, int nextCol,
                          const MatrixXd &H) {
  std::ofstream out(path);
  if (!out) {
    return false;
  }
  out << "eon_hess_ckpt " << size << " " << nextCol << "\n";
  out.precision(17);
  for (int i = 0; i < size; ++i) {
    for (int j = 0; j < size; ++j) {
      out << H(i, j) << (j + 1 == size ? '\n' : ' ');
    }
  }
  return static_cast<bool>(out);
}

} // namespace

Hessian::Hessian(const Parameters &params, Matter *matter)
    : matter{matter},
      parameters{params} {
  hessian.resize(0, 0);
  freqs.resize(0);
}

MatrixXd Hessian::getHessian(Matter *matterIn, const VectorXi &atomsIn) {
  if ((matter != matterIn) || (atoms.size() != atomsIn.size()) ||
      (atoms != atomsIn) || (hessian.rows() == 0)) {
    hessian.resize(0, 0);
    matter = matterIn;
    atoms = atomsIn;

    if (!calculate()) {
      hessian.resize(0, 0);
    }
  }
  return hessian;
}

VectorXd Hessian::getFreqs(Matter *matterIn, const VectorXi &atomsIn) {
  if ((matter != matterIn) || (atoms.size() != atomsIn.size()) ||
      (atoms != atomsIn) || (hessian.rows() == 0)) {
    hessian.resize(0, 0);
    matter = matterIn;
    atoms = atomsIn;

    if (!calculate()) {
      freqs.resize(0);
    }
  }
  return freqs;
}

bool Hessian::calculate() {
  int nAtoms = matter->numberOfAtoms();

  int size = static_cast<int>(atoms.rows()) * 3;
  QUILL_LOG_DEBUG(log, "[Hessian] Hessian size: {}\n", size);
  if (size == 0) {
    return false;
  }

  // Mobile-atom polarity: indices in `atoms` are FD-displaced DOF owners.
  for (int a = 0; a < atoms.rows(); ++a) {
    const long idx = atoms(a);
    if (idx < 0 || idx >= nAtoms) {
      QUILL_LOG_ERROR(log,
                      "[Hessian] atom index {} out of range [0, {}) at list "
                      "entry {}; aborting FD Hessian",
                      idx, nAtoms, a);
      return false;
    }
  }

  Matter matterTemp(*matter);
  double dr = parameters.main_options.finiteDifference;
  if (!(dr > 0.0) || !std::isfinite(dr)) {
    QUILL_LOG_ERROR(log, "[Hessian] invalid finiteDifference dr={}\n", dr);
    return false;
  }

  const bool useCentral =
      isCentralScheme(parameters.hessian_options.fd_scheme);
  const std::string &ckptPath = parameters.hessian_options.checkpoint_path;
  const bool wantResume = parameters.hessian_options.resume && !ckptPath.empty();

  AtomMatrix pos = matter->getPositions();
  AtomMatrix posDisplace(nAtoms, 3);
  AtomMatrix posTemp(nAtoms, 3);
  AtomMatrix force0(nAtoms, 3);
  AtomMatrix forcePlus(nAtoms, 3);
  AtomMatrix forceMinus(nAtoms, 3);

  hessian.resize(size, size);
  hessian.setZero();

  int startCol = 0;
  if (wantResume && loadColumnCheckpoint(ckptPath, size, startCol, hessian)) {
    QUILL_LOG_DEBUG(log, "[Hessian] resume from column {} / {}\n", startCol,
                    size);
  } else {
    startCol = 0;
    hessian.setZero();
  }

  force0 = matterTemp.getForces();
  if (!force0.allFinite()) {
    QUILL_LOG_ERROR(log,
                    "[Hessian] non-finite forces at undisplaced geometry; "
                    "aborting FD Hessian");
    return false;
  }

  for (int i = startCol; i < size; i++) {
    posDisplace.setZero();
    posDisplace(atoms(i / 3), i % 3) = dr;

    posTemp = pos + posDisplace;
    matterTemp.setPositions(posTemp);
    forcePlus = matterTemp.getForces();
    if (!forcePlus.allFinite()) {
      QUILL_LOG_ERROR(log,
                      "[Hessian] non-finite forces for FD column {} (+); "
                      "aborting FD Hessian",
                      i);
      return false;
    }

    if (useCentral) {
      posTemp = pos - posDisplace;
      matterTemp.setPositions(posTemp);
      forceMinus = matterTemp.getForces();
      if (!forceMinus.allFinite()) {
        QUILL_LOG_ERROR(log,
                        "[Hessian] non-finite forces for FD column {} (-); "
                        "aborting FD Hessian",
                        i);
        return false;
      }
      // Central: H_ij ≈ -(F+(xj) - F-(xj)) / (2 dr), mass-weighted
      for (int j = 0; j < size; j++) {
        const double dF =
            forcePlus(atoms(j / 3), j % 3) - forceMinus(atoms(j / 3), j % 3);
        hessian(i, j) = -dF / (2.0 * dr);
        const double effMass = std::sqrt(matter->getMass(atoms(j / 3)) *
                                         matter->getMass(atoms(i / 3)));
        hessian(i, j) = eonc::safemath::safe_div(hessian(i, j), effMass, 0.0);
      }
    } else {
      // One-sided (forward): H_ij ≈ -(F+(xj) - F0(xj)) / dr  [default; cheaper]
      for (int j = 0; j < size; j++) {
        const double dF =
            forcePlus(atoms(j / 3), j % 3) - force0(atoms(j / 3), j % 3);
        hessian(i, j) = -dF / dr;
        const double effMass = std::sqrt(matter->getMass(atoms(j / 3)) *
                                         matter->getMass(atoms(i / 3)));
        hessian(i, j) = eonc::safemath::safe_div(hessian(i, j), effMass, 0.0);
      }
    }

    if (!ckptPath.empty()) {
      // next column to compute after a clean interrupt
      saveColumnCheckpoint(ckptPath, size, i + 1, hessian);
    }
  }

  // Symmetrize (FD noise breaks H=H^T; required for vib analysis)
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < i; j++) {
      hessian(i, j) = (hessian(i, j) + hessian(j, i)) / 2;
      hessian(j, i) = hessian(i, j);
    }
  }

  if (!hessian.allFinite()) {
    QUILL_LOG_ERROR(log,
                    "[Hessian] non-finite entries after FD assembly; "
                    "aborting eigen solve");
    return false;
  }

  if (!parameters.main_options.quiet) {
    QUILL_LOG_DEBUG(log, "[Hessian] writing hessian\n");
    std::ofstream hessfile;
    hessfile.open("hessian.dat");
    hessfile << hessian;
    hessfile.close();
  }

  // Completed run: remove checkpoint so a later job does not resume stale cols
  if (!ckptPath.empty()) {
    std::remove(ckptPath.c_str());
  }

  double t0, t1;
  eonc::helpers::getTime(&t0, nullptr, nullptr);
  QUILL_LOG_DEBUG(log, "[Hessian] calculating eigen values of the hessian\n");
  // ColMajor copy for SelfAdjointEigenSolver (eOn MatrixXd is RowMajor)
  using ColMajorXd =
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
  ColMajorXd hessianCol = hessian;
  Eigen::SelfAdjointEigenSolver<ColMajorXd> es(hessianCol,
                                               Eigen::EigenvaluesOnly);
  eonc::helpers::getTime(&t1, nullptr, nullptr);
  QUILL_LOG_DEBUG(log, "[Hessian] eigenvalue problem took {:.4e} seconds\n",
                  t1 - t0);
  if (es.info() != Eigen::Success) {
    QUILL_LOG_ERROR(log,
                    "[Hessian] SelfAdjointEigenSolver failed (info={}); "
                    "aborting",
                    static_cast<int>(es.info()));
    return false;
  }
  freqs = es.eigenvalues();
  if (!freqs.allFinite()) {
    QUILL_LOG_ERROR(log, "[Hessian] non-finite eigenvalues; aborting");
    return false;
  }

  return true;
}

VectorXd Hessian::removeZeroFreqs(const VectorXd &freqs) {
  QUILL_LOG_DEBUG(log, "[Hessian] removing zero frequency modes");
  int size = freqs.size();
  if (size != 3 * matter->numberOfAtoms()) {
    return freqs;
  }
  VectorXd newfreqs;
  newfreqs.resize(size);
  int nremoved = 0;
  for (int i = 0; i < size; i++) {
    if (std::abs(freqs(i)) > parameters.hessian_options.zero_freq_value) {
      newfreqs(i - nremoved) = freqs(i);
    } else {
      nremoved++;
    }
  }

  if (nremoved != 6) {
    QUILL_LOG_ERROR(
        log, "[Hessian] [error] Found {} trivial eigenmodes instead of 6",
        nremoved);
  }
  return newfreqs.head(size - nremoved);
}
