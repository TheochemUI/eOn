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
#include "eon/Hessian.h"
#include "eon/EonLogger.h"
#include "eon/HelperFunctions.h"
#include "eon/SafeMath.h"
#include "eon/VesinNeighbors.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace eonc {

namespace {

// phva_atoms entries are *mobile / displaced* atoms for FD (hybrid/PHVA-class
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
      hessian.resize(0, 0);
    }
  }
  return freqs;
}

namespace {

struct MobileColoring {
  std::vector<int> color;
  // closed[ia] = local mobile indices in the closed cutoff neighborhood of ia
  std::vector<std::vector<int>> closed;
};

MobileColoring buildMobileColoring(const Matter &matter, const VectorXi &atoms,
                                   double cutoff) {
  MobileColoring out;
  const int nAtoms = static_cast<int>(matter.numberOfAtoms());
  const int nMobile = static_cast<int>(atoms.rows());
  if (nAtoms <= 0 || nMobile <= 0 || !(cutoff > 0.0) ||
      !std::isfinite(cutoff)) {
    return out;
  }
  const Matrix3d cell = matter.getCell();
  if (!(std::abs(cell.determinant()) > 1e-18)) {
    return out;
  }
  for (int a = 0; a < nMobile; ++a) {
    const long idx = atoms(a);
    if (idx < 0 || idx >= nAtoms) {
      return {};
    }
  }

  const AtomMatrix &pos = matter.getPositions();
  VesinNeighbors nl;
  VesinNeighbors::Options opt;
  opt.cutoff = cutoff;
  opt.full = true;
  opt.return_distances = false;
  opt.return_vectors = false;
  opt.periodic = {matter.getPeriodic(), matter.getPeriodic(),
                  matter.getPeriodic()};
  try {
    nl.compute(pos.data(), static_cast<std::size_t>(nAtoms), cell.data(), opt);
  } catch (const std::exception &) {
    return {};
  }

  std::vector<std::vector<int>> nbs(static_cast<std::size_t>(nAtoms));
  for (std::size_t p = 0; p < nl.size(); ++p) {
    const int i = static_cast<int>(nl.i(p));
    const int j = static_cast<int>(nl.j(p));
    if (i == j || i < 0 || j < 0 || i >= nAtoms || j >= nAtoms) {
      continue;
    }
    nbs[static_cast<std::size_t>(i)].push_back(j);
  }
  for (auto &row : nbs) {
    std::sort(row.begin(), row.end());
    row.erase(std::unique(row.begin(), row.end()), row.end());
  }

  std::vector<int> local(static_cast<std::size_t>(nAtoms), -1);
  for (int a = 0; a < nMobile; ++a) {
    local[static_cast<std::size_t>(atoms(a))] = a;
  }

  out.closed.assign(static_cast<std::size_t>(nMobile), {});
  for (int ia = 0; ia < nMobile; ++ia) {
    auto &nbhd = out.closed[static_cast<std::size_t>(ia)];
    nbhd.push_back(ia);
    const int g = static_cast<int>(atoms(ia));
    for (int nb : nbs[static_cast<std::size_t>(g)]) {
      const int lb = local[static_cast<std::size_t>(nb)];
      if (lb >= 0) {
        nbhd.push_back(lb);
      }
    }
    std::sort(nbhd.begin(), nbhd.end());
    nbhd.erase(std::unique(nbhd.begin(), nbhd.end()), nbhd.end());
  }

  // Square of the cutoff graph on the mobile set: atoms whose closed
  // neighborhoods intersect cannot move in the same finite-difference.
  std::vector<std::vector<int>> touchers(static_cast<std::size_t>(nMobile));
  for (int ia = 0; ia < nMobile; ++ia) {
    for (int k : out.closed[static_cast<std::size_t>(ia)]) {
      touchers[static_cast<std::size_t>(k)].push_back(ia);
    }
  }
  std::vector<std::vector<int>> adj(static_cast<std::size_t>(nMobile));
  for (int k = 0; k < nMobile; ++k) {
    const auto &t = touchers[static_cast<std::size_t>(k)];
    for (size_t a = 0; a < t.size(); ++a) {
      for (size_t b = a + 1; b < t.size(); ++b) {
        adj[static_cast<std::size_t>(t[a])].push_back(t[b]);
        adj[static_cast<std::size_t>(t[b])].push_back(t[a]);
      }
    }
  }
  for (auto &row : adj) {
    std::sort(row.begin(), row.end());
    row.erase(std::unique(row.begin(), row.end()), row.end());
  }
  out.color = greedyColorCutoffGraph(adj);
  return out;
}

void writeMassWeighted(MatrixXd &hessian, const Matter &matter,
                       const VectorXi &atoms, int col, int atomI,
                       const AtomMatrix &forceA, const AtomMatrix &forceB,
                       double denom, const std::vector<int> &owner, int ia) {
  const int size = static_cast<int>(atoms.rows()) * 3;
  const double massI = matter.getMass(atomI);
  for (int j = 0; j < size; ++j) {
    const long atomJ = atoms(j / 3);
    double dF = 0.0;
    if (owner[static_cast<std::size_t>(j / 3)] == ia) {
      dF = forceA(atomJ, j % 3) - forceB(atomJ, j % 3);
    }
    hessian(col, j) = -dF / denom;
    const double effMass = std::sqrt(matter.getMass(atomJ) * massI);
    hessian(col, j) = eonc::safemath::safe_div(hessian(col, j), effMass, 0.0);
  }
}

} // namespace

std::vector<int>
greedyColorCutoffGraph(const std::vector<std::vector<int>> &adj) {
  const int n = static_cast<int>(adj.size());
  std::vector<int> color(static_cast<std::size_t>(n), -1);
  std::vector<int> used(static_cast<std::size_t>(std::max(n, 0)), 0);
  int epoch = 0;
  for (int v = 0; v < n; ++v) {
    ++epoch;
    for (int u : adj[static_cast<std::size_t>(v)]) {
      if (u < 0 || u >= n || u == v) {
        continue;
      }
      const int cu = color[static_cast<std::size_t>(u)];
      if (cu >= 0 && cu < n) {
        used[static_cast<std::size_t>(cu)] = epoch;
      }
    }
    int c = 0;
    while (c < n && used[static_cast<std::size_t>(c)] == epoch) {
      ++c;
    }
    color[static_cast<std::size_t>(v)] = c;
  }
  return color;
}

std::vector<int> colorMobileCutoffGraph(const Matter &matter,
                                        const VectorXi &atoms, double cutoff) {
  return buildMobileColoring(matter, atoms, cutoff).color;
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

  double dr = parameters.main_options().finiteDifference;
  if (!(dr > 0.0) || !std::isfinite(dr)) {
    QUILL_LOG_ERROR(log, "[Hessian] invalid finiteDifference dr={}\n", dr);
    return false;
  }

  const bool useCentral =
      isCentralScheme(parameters.hessian_options().fd_scheme);
  const std::string &ckptPath = parameters.hessian_options().checkpoint_path;

  hessian.resize(size, size);
  hessian.setZero();

  // Net-force removal adds the same shift to every atom, so columns are
  // no longer confined to the cutoff neighborhood. A column checkpoint is
  // also stored one coordinate at a time.
  bool anyFixed = false;
  for (int i = 0; i < nAtoms; ++i) {
    if (matter->getFixed(i)) {
      anyFixed = true;
      break;
    }
  }
  const bool netCoupled =
      parameters.main_options().removeNetForce && nAtoms > 1 && !anyFixed;
  double cutoff = 0.0;
  if (matter->getPotential()) {
    cutoff = matter->getPotential()->finiteCutoff();
  }
  if (!netCoupled && ckptPath.empty() && cutoff > 0.0 &&
      std::isfinite(cutoff)) {
    if (calculateColored(cutoff, dr, useCentral)) {
      return true;
    }
    hessian.setZero();
  }
  return calculateSerial(dr, useCentral);
}

bool Hessian::calculateColored(double cutoff, double dr, bool useCentral) {
  const int nAtoms = static_cast<int>(matter->numberOfAtoms());
  const int nMobile = static_cast<int>(atoms.rows());
  const int size = nMobile * 3;
  const MobileColoring coloring = buildMobileColoring(*matter, atoms, cutoff);
  if (static_cast<int>(coloring.color.size()) != nMobile ||
      static_cast<int>(coloring.closed.size()) != nMobile) {
    return false;
  }

  int nColors = 0;
  for (int c : coloring.color) {
    if (c < 0) {
      return false;
    }
    nColors = std::max(nColors, c + 1);
  }
  if (nColors <= 0) {
    return false;
  }
  QUILL_LOG_DEBUG(log,
                  "[Hessian] cutoff coloring: {} colors for {} mobile atoms\n",
                  nColors, nMobile);

  std::vector<std::vector<int>> members(static_cast<std::size_t>(nColors));
  for (int ia = 0; ia < nMobile; ++ia) {
    members[static_cast<std::size_t>(
                coloring.color[static_cast<std::size_t>(ia)])]
        .push_back(ia);
  }
  std::vector<std::vector<int>> owner(
      static_cast<std::size_t>(nColors),
      std::vector<int>(static_cast<std::size_t>(nMobile), -1));
  for (int c = 0; c < nColors; ++c) {
    for (int ia : members[static_cast<std::size_t>(c)]) {
      for (int k : coloring.closed[static_cast<std::size_t>(ia)]) {
        int &slot =
            owner[static_cast<std::size_t>(c)][static_cast<std::size_t>(k)];
        if (slot >= 0 && slot != ia) {
          return false;
        }
        slot = ia;
      }
    }
  }

  Matter matterTemp(*matter);
  const AtomMatrix pos = matter->getPositions();
  AtomMatrix posDisplace(nAtoms, 3);
  AtomMatrix force0 = matterTemp.getForces();
  if (!force0.allFinite()) {
    QUILL_LOG_ERROR(log, "[Hessian] non-finite forces at undisplaced geometry; "
                         "aborting FD Hessian");
    return false;
  }

  for (int dir = 0; dir < 3; ++dir) {
    for (int c = 0; c < nColors; ++c) {
      const auto &group = members[static_cast<std::size_t>(c)];
      if (group.empty()) {
        continue;
      }
      posDisplace.setZero();
      for (int ia : group) {
        posDisplace(atoms(ia), dir) = dr;
      }
      matterTemp.setPositions(pos + posDisplace);
      AtomMatrix forcePlus = matterTemp.getForces();
      if (!forcePlus.allFinite()) {
        QUILL_LOG_ERROR(log,
                        "[Hessian] non-finite forces for color {} dir {} (+); "
                        "aborting FD Hessian",
                        c, dir);
        return false;
      }
      AtomMatrix forceMinus;
      if (useCentral) {
        matterTemp.setPositions(pos - posDisplace);
        forceMinus = matterTemp.getForces();
        if (!forceMinus.allFinite()) {
          QUILL_LOG_ERROR(log,
                          "[Hessian] non-finite forces for color {} dir {} "
                          "(-); aborting FD Hessian",
                          c, dir);
          return false;
        }
      }
      const double denom = useCentral ? (2.0 * dr) : dr;
      const AtomMatrix &forceB = useCentral ? forceMinus : force0;
      for (int ia : group) {
        const int col = ia * 3 + dir;
        writeMassWeighted(hessian, *matter, atoms, col,
                          static_cast<int>(atoms(ia)), forcePlus, forceB, denom,
                          owner[static_cast<std::size_t>(c)], ia);
      }
    }
  }
  return finalizeHessian(size);
}

bool Hessian::calculateSerial(double dr, bool useCentral) {
  const int nAtoms = static_cast<int>(matter->numberOfAtoms());
  const int size = static_cast<int>(atoms.rows()) * 3;
  const std::string &ckptPath = parameters.hessian_options().checkpoint_path;
  const bool wantResume =
      parameters.hessian_options().resume && !ckptPath.empty();

  AtomMatrix pos = matter->getPositions();
  AtomMatrix posDisplace(nAtoms, 3);
  AtomMatrix posTemp(nAtoms, 3);
  AtomMatrix forcePlus(nAtoms, 3);
  AtomMatrix forceMinus(nAtoms, 3);

  int startCol = 0;
  if (wantResume && loadColumnCheckpoint(ckptPath, size, startCol, hessian)) {
    QUILL_LOG_DEBUG(log, "[Hessian] resume from column {} / {}\n", startCol,
                    size);
  } else {
    startCol = 0;
    hessian.setZero();
  }

  Matter matterTemp(*matter);
  AtomMatrix force0 = matterTemp.getForces();
  if (!force0.allFinite()) {
    QUILL_LOG_ERROR(log, "[Hessian] non-finite forces at undisplaced geometry; "
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
      // One-sided (forward): H_ij ≈ -(F+(xj) - F0(xj)) / dr
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
  return finalizeHessian(size);
}

bool Hessian::finalizeHessian(int size) {
  const std::string &ckptPath = parameters.hessian_options().checkpoint_path;

  // Symmetrize (FD noise breaks H=H^T; required for vib analysis)
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < i; j++) {
      hessian(i, j) = (hessian(i, j) + hessian(j, i)) / 2;
      hessian(j, i) = hessian(i, j);
    }
  }

  if (!hessian.allFinite()) {
    QUILL_LOG_ERROR(log, "[Hessian] non-finite entries after FD assembly; "
                         "aborting eigen solve");
    return false;
  }

  if (!parameters.main_options().quiet) {
    QUILL_LOG_DEBUG(log, "[Hessian] writing hessian\n");
  }
  {
    std::ofstream hessfile("hessian.dat");
    if (!hessfile) {
      QUILL_LOG_ERROR(log, "[Hessian] failed to open hessian.dat");
      return false;
    }
    hessfile << hessian;
    hessfile.close();
    if (!hessfile) {
      QUILL_LOG_ERROR(log, "[Hessian] failed to write hessian.dat");
      return false;
    }
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
    if (std::abs(freqs(i)) > parameters.hessian_options().zero_freq_value) {
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

} // namespace eonc
