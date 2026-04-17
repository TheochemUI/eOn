/*
 * This file is part of eOn.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 *
 * Copyright (c) 2010--present, eOn Development Team
 * All rights reserved.
 *
 * Repo:
 * https://github.com/TheochemUI/eOn
 */
#include "IRACompare.h"
#include "Eigen.h"
#include <cstdlib>

#ifdef WITH_IRA
extern "C" {
#include "iralib_interf.h"
}
#include "libs/IRA/IRAResource.h" // Include IRA resource after IRA interfaces are defined
#endif

namespace eonc {

IRACompare::MatchResult IRACompare::match(const Matter &m1, const Matter &m2,
                                          double distThreshold) {
  MatchResult result;
#ifdef WITH_IRA
  auto &res = get_ira_resource();

  // Lock the library for the duration of this specific comparison
  std::lock_guard<std::mutex> lock(res.library_mutex);

  try {
    res.require_loaded();
  } catch (const std::exception &e) {
    result.error = -1;
    return result;
  }

  const int nat1 = m1.numberOfAtoms();
  const int nat2 = m2.numberOfAtoms();

  // Prepare type arrays
  std::vector<int> typ1(nat1), typ2(nat2);
  auto nrs1 = m1.getAtomicNrs();
  auto nrs2 = m2.getAtomicNrs();
  for (int i = 0; i < nat1; i++)
    typ1[i] = nrs1[i];
  for (int i = 0; i < nat2; i++)
    typ2[i] = nrs2[i];

  // Prepare coordinates using the new Eigen helper - use direct data access
  const AtomMatrix &pos1 = m1.getPositions();
  const AtomMatrix &pos2 = m2.getPositions();

  // Use Eigen::Map to reinterpret the row-major data as column-major for
  // Fortran
  Eigen::Map<const AtomMatrixF> coords1_map(pos1.data(), 3, nat1);
  Eigen::Map<const AtomMatrixF> coords2_map(pos2.data(), 3, nat2);

  // Candidate arrays: -1 means "use geometric center as origin" (good
  // initial guess for translation with equal-size structures)
  std::vector<int> cand1(nat1, 0), cand2(nat2, 0);
  if (nat1 == nat2) {
    cand1[0] = -1;
    cand2[0] = -1;
  } else {
    cand1[0] = 1;
    for (int i = 0; i < nat2; i++)
      cand2[i] = i + 1;
  }

  // Pre-allocate output buffers (libira_match fills these)
  std::vector<double> rmat_buf(9);
  std::vector<double> tr_buf(3);
  std::vector<int> perm_buf(nat2);
  double hd = 0.0;
  int ierr = 0;

  // libira_match expects double** for rotation and translation
  // (it may reallocate internally)
  double *rmat_ptr = rmat_buf.data();
  double *tr_ptr = tr_buf.data();
  int *perm_ptr = perm_buf.data();

  res.get_match_fn()(nat1, typ1.data(), coords1_map.data(), cand1.data(), nat2,
                     typ2.data(), coords2_map.data(), cand2.data(),
                     distThreshold, &rmat_ptr, &tr_ptr, &perm_ptr, &hd, &ierr);

  result.error = ierr;
  if (ierr == 0) {
    result.permutation.assign(perm_ptr, perm_ptr + nat2);
    result.hausdorffDistance = hd;
    // rmat is 3x3 flat array (column-major from Fortran)
    // Use Eigen maps for efficient conversion
    Eigen::Map<const Matrix3d> rot_map(rmat_ptr);
    result.rotation =
        rot_map.transpose(); // Transpose to get proper row-major layout
    result.translation = Eigen::Map<const Vector3d>(tr_ptr);
  }

  // If libira reallocated (pointers changed), free the new buffers.
  // libira allocates replacement output buffers via c_malloc on the Fortran
  // side (see libira's lib_match wrapper), so std::free is the matching
  // deallocator. Never free the original stack/vector-backed buffers.
  if (rmat_ptr != rmat_buf.data())
    std::free(rmat_ptr);
  if (tr_ptr != tr_buf.data())
    std::free(tr_ptr);
  if (perm_ptr != perm_buf.data())
    std::free(perm_ptr);
#else
  result.error = -1;
#endif
  return result;
}

IRACompare::MatchResult IRACompare::matchPBC(const Matter &m1, const Matter &m2,
                                             double distThreshold) {
  MatchResult result;
#ifdef WITH_IRA
  auto &res = get_ira_resource();

  // Lock the library for the duration of this specific comparison
  std::lock_guard<std::mutex> lock(res.library_mutex);

  try {
    res.require_loaded();
  } catch (const std::exception &e) {
    result.error = -1;
    return result;
  }

  const int nat1 = m1.numberOfAtoms();
  const int nat2 = m2.numberOfAtoms();

  std::vector<int> typ1(nat1), typ2(nat2);
  auto nrs1 = m1.getAtomicNrs();
  auto nrs2 = m2.getAtomicNrs();
  for (int i = 0; i < nat1; i++)
    typ1[i] = nrs1[i];
  for (int i = 0; i < nat2; i++)
    typ2[i] = nrs2[i];

  // Prepare coordinates using direct Eigen data access
  const AtomMatrix &pos1 = m1.getPositions();
  const AtomMatrix &pos2 = m2.getPositions();

  Eigen::Map<const AtomMatrixF> coords1_map(pos1.data(), 3, nat1);
  Eigen::Map<const AtomMatrixF> coords2_map(pos2.data(), 3, nat2);

  // Lattice vectors (3x3 column-major)
  Matrix3d cell = m2.getCell();
  double lat[9];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      lat[j * 3 + i] = cell(i, j);

  // Caller must pre-allocate output arrays; Fortran receives pointers only
  std::vector<int> found_buf(nat1);
  std::vector<double> dists_buf(nat1);
  int *found_ptr = found_buf.data();
  double *dists_ptr = dists_buf.data();

  // NOTE: this calls cshda_pbc only (assignment, no rotation/SVD matching).
  // It finds the best atom assignment under PBC but does not compute the
  // optimal rotation or translation.
  res.get_cshda_pbc_fn()(nat1, typ1.data(), coords1_map.data(), nat2,
                         typ2.data(), coords2_map.data(), lat, distThreshold,
                         &found_ptr, &dists_ptr);

  result.permutation.assign(found_ptr, found_ptr + nat1);
  result.hausdorffDistance = 0.0;
  for (int i = 0; i < nat1; i++) {
    result.hausdorffDistance = std::max(result.hausdorffDistance, dists_ptr[i]);
  }
  result.rotation = Eigen::Matrix3d::Identity();
  result.translation = Eigen::Vector3d::Zero();
  result.error = 0;
#else
  result.error = -1;
#endif
  return result;
}

IRACompare::SymmetryResult
IRACompare::findSymmetry(const Matter &m, double threshold, bool prescreenIh) {
  SymmetryResult result;
#ifdef WITH_IRA
  auto &res = get_ira_resource();

  // Lock the library for the duration of this specific symmetry analysis
  std::lock_guard<std::mutex> lock(res.library_mutex);

  try {
    res.require_loaded();
  } catch (const std::exception &e) {
    result.error = -1;
    return result;
  }

  const int nat = m.numberOfAtoms();
  std::vector<int> typ(nat);
  auto nrs = m.getAtomicNrs();
  for (int i = 0; i < nat; i++)
    typ[i] = nrs[i];

  // Prepare coordinates using direct Eigen data access
  const AtomMatrix &pos = m.getPositions();
  Eigen::Map<const AtomMatrixF> coords_map(pos.data(), 3, nat);

  // libira_compute_all requires pre-allocated output arrays of size nmax
  int nmax = res.get_get_nmax_fn()();

  int n_mat = 0;
  std::vector<double> mat_buf(9 * nmax);
  std::vector<int> perm_buf(nat * nmax);
  std::vector<char> op_buf(nmax + 1);
  std::vector<int> n_buf(nmax);
  std::vector<int> p_buf(nmax);
  std::vector<double> ax_buf(3 * nmax);
  std::vector<double> angle_buf(nmax);
  std::vector<double> dH_buf(nmax);
  std::vector<char> pg_buf(11);
  int n_prin_ax = 0;
  std::vector<double> prin_ax_buf(3 * nmax);
  int cerr = 0;

  double *mat_data = mat_buf.data();
  int *perm_data = perm_buf.data();
  char *op_data = op_buf.data();
  int *n_data = n_buf.data();
  int *p_data = p_buf.data();
  double *ax_data = ax_buf.data();
  double *angle_data = angle_buf.data();
  double *dH_data = dH_buf.data();
  char *pg = pg_buf.data();
  double *prin_ax = prin_ax_buf.data();

  res.get_compute_all_fn()(nat, typ.data(), coords_map.data(), threshold,
                           prescreenIh ? 1 : 0, &n_mat, &mat_data, &perm_data,
                           &op_data, &n_data, &p_data, &ax_data, &angle_data,
                           &dH_data, &pg, &n_prin_ax, &prin_ax, &cerr);

  result.error = cerr;
  if (cerr == 0) {
    result.nOperations = n_mat;
    if (pg) {
      result.pointGroup = std::string(pg);
    }

    result.operations.resize(n_mat);
    for (int k = 0; k < n_mat; k++) {
      Eigen::Map<const Matrix3d> op_map(&mat_data[k * 9]);
      result.operations[k] =
          op_map.transpose(); // Transpose to get proper row-major layout
    }

    if (angle_data) {
      result.angles.assign(angle_data, angle_data + n_mat);
    }
    if (ax_data) {
      result.axes.resize(n_mat);
      for (int k = 0; k < n_mat; k++) {
        result.axes[k] = Eigen::Map<const Vector3d>(&ax_data[k * 3]);
      }
    }
  }

  // Only free if Fortran reallocated (pointer changed from our pre-allocated
  // buffer). If pointer is unchanged, the vector destructor handles cleanup.
  // Reallocated buffers come from libira's c_malloc side, so std::free is the
  // matching deallocator; never free the vector-backed originals.
  if (mat_data != mat_buf.data())
    std::free(mat_data);
  if (perm_data != perm_buf.data())
    std::free(perm_data);
  if (op_data != op_buf.data())
    std::free(op_data);
  if (n_data != n_buf.data())
    std::free(n_data);
  if (p_data != p_buf.data())
    std::free(p_data);
  if (ax_data != ax_buf.data())
    std::free(ax_data);
  if (angle_data != angle_buf.data())
    std::free(angle_data);
  if (dH_data != dH_buf.data())
    std::free(dH_data);
  if (pg != pg_buf.data())
    std::free(pg);
  if (prin_ax != prin_ax_buf.data())
    std::free(prin_ax);
#else
  result.error = -1;
#endif
  return result;
}

} // namespace eonc
