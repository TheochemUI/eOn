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

Hessian::Hessian(const Parameters &params, Matter *matter)
    : matter{matter},
      parameters{params} {
  hessian.resize(0, 0);
  freqs.resize(0);
  /* Logger initialized via class member */
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

  // Determine the Hessian size
  int size = 0;
  size = static_cast<int>(atoms.rows()) * 3;
  QUILL_LOG_DEBUG(log, "[Hessian] Hessian size: {}\n", size);
  if (size == 0) {
    return false;
  }

  // Out-of-range atom indices used to OOB-write AtomMatrix rows and segfault
  // inside potential force evaluation / Eigen (seen with VTST partial Hessians
  // and portable eonclient builds). Fail cleanly instead.
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

  // Build the hessian
  Matter matterTemp(*matter);
  double dr = parameters.main_options.finiteDifference;
  if (!(dr > 0.0) || !std::isfinite(dr)) {
    QUILL_LOG_ERROR(log, "[Hessian] invalid finiteDifference dr={}\n", dr);
    return false;
  }

  AtomMatrix pos = matter->getPositions();
  AtomMatrix posDisplace(nAtoms, 3);
  AtomMatrix posTemp(nAtoms, 3);
  AtomMatrix force1(nAtoms, 3);
  AtomMatrix force2(nAtoms, 3);

  //    Matrix <double, Eigen::Dynamic, Eigen::Dynamic> hessian(size, size);
  hessian.resize(size, size);

  force1 = matterTemp.getForces();
  if (!force1.allFinite()) {
    QUILL_LOG_ERROR(log,
                    "[Hessian] non-finite forces at undisplaced geometry; "
                    "aborting FD Hessian");
    return false;
  }
  for (int i = 0; i < size; i++) {
    posDisplace.setZero();

    // Displacing one coordinate
    posDisplace(atoms(i / 3), i % 3) = dr;

    posTemp = pos + posDisplace;
    matterTemp.setPositions(posTemp);
    force2 = matterTemp.getForces();
    if (!force2.allFinite()) {
      QUILL_LOG_ERROR(log,
                      "[Hessian] non-finite forces for FD column {}; "
                      "aborting FD Hessian",
                      i);
      return false;
    }

    // To use central difference estimate of the hessian uncomment following
    // (and divide by 2*dr) in the following. This does use an additional 'size'
    // forcecalls and will generally not lead to very different results. In most
    // cases, the additional accuracy is not worth the computation time.

    /*
    posTemp = pos - posDisplace;
    matterTemp.setPositions(posTemp);
    force1 = matterTemp.getForces();
    */

    for (int j = 0; j < size; j++) {
      hessian(i, j) =
          -(force2(atoms(j / 3), j % 3) - force1(atoms(j / 3), j % 3)) / dr;
      double effMass = std::sqrt(matter->getMass(atoms(j / 3)) *
                                 matter->getMass(atoms(i / 3)));
      hessian(i, j) = eonc::safemath::safe_div(hessian(i, j), effMass, 0.0);
    }
  }

  // Force hessian to be symmetric

  // hessian = (hessian + hessian.transpose())/2;
  // cannot be used, messes up the lower trianguler
  // transpose does not seem to be a hardcopy, rather just an index manipulation

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

  double t0, t1;
  eonc::helpers::getTime(&t0, nullptr, nullptr);
  QUILL_LOG_DEBUG(log, "[Hessian] calculating eigen values of the hessian\n");
  // eOn's MatrixXd is RowMajor (atom-block layout). Eigen's
  // SelfAdjointEigenSolver is only well-defined for ColMajor storage; running
  // it on RowMajor has caused segfaults in HessianJob (VTST partial blocks,
  // size ~42) on some toolchains. Copy to ColMajor for the eigen solve only.
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

// If we are checking for rotation, then the system has no frozen atoms and
// can rotate and translate. This gives effectively zero eigenvalues. We
// need to remove them from the prefactor calculation.
// the condition requires that every atom moves. Otherwise, we don't
// get the 6 rotational and translational modes.
// XXX: what happens if the entire particle rotates about one atom or a line of
// atoms?

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
