#pragma once
/**
 * Helpers: eOn row-major Eigen matrices <-> numpy ndarrays (nanobind).
 *
 * Prefer *view* helpers for Matter-owned storage (zero-copy, keep parent
 * alive via rv_policy::reference_internal). Use *copy* helpers when the
 * Eigen temporary must outlive the Python buffer.
 */
#include "eon/Eigen.h"

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include <cstring>
#include <stdexcept>
#include <string>

namespace eonc::pybind {
namespace nb = nanobind;

using NpF64 = nb::ndarray<nb::numpy, double, nb::c_contig, nb::device::cpu>;
using NpI32 = nb::ndarray<nb::numpy, int32_t, nb::c_contig, nb::device::cpu>;
using NpI64 = nb::ndarray<nb::numpy, int64_t, nb::c_contig, nb::device::cpu>;

/// Zero-copy (n, 3) float64 view of external/row-major AtomMatrix storage.
inline nb::ndarray<nb::numpy, double, nb::c_contig, nb::device::cpu>
view_n3(double *data, long n) {
  if (n < 0) {
    throw std::invalid_argument("view_n3: negative n");
  }
  return nb::ndarray<nb::numpy, double, nb::c_contig, nb::device::cpu>(
      data, {static_cast<size_t>(n), size_t{3}});
}

inline nb::ndarray<nb::numpy, double, nb::c_contig, nb::device::cpu>
view_n3(const double *data, long n) {
  return view_n3(const_cast<double *>(data), n);
}

inline nb::ndarray<nb::numpy, double, nb::c_contig, nb::device::cpu>
view_33(double *data) {
  return nb::ndarray<nb::numpy, double, nb::c_contig, nb::device::cpu>(
      data, {size_t{3}, size_t{3}});
}

inline nb::ndarray<nb::numpy, double, nb::c_contig, nb::device::cpu>
view_n(double *data, long n) {
  return nb::ndarray<nb::numpy, double, nb::c_contig, nb::device::cpu>(
      data, {static_cast<size_t>(n)});
}

inline nb::ndarray<nb::numpy, int32_t, nb::c_contig, nb::device::cpu>
view_n_i32(int *data, long n) {
  static_assert(sizeof(int) == 4, "view_n_i32 assumes 32-bit int");
  return nb::ndarray<nb::numpy, int32_t, nb::c_contig, nb::device::cpu>(
      reinterpret_cast<int32_t *>(data), {static_cast<size_t>(n)});
}

/// Require C-contiguous (n,3) float64; return data pointer (no copy).
inline const double *require_n3(const NpF64 &arr, long expected_n = -1) {
  if (arr.ndim() != 2 || arr.shape(1) != 3) {
    throw std::invalid_argument(
        "expected float64 array of shape (n, 3), got ndim=" +
        std::to_string(arr.ndim()));
  }
  if (expected_n >= 0 && static_cast<long>(arr.shape(0)) != expected_n) {
    throw std::invalid_argument("expected n=" + std::to_string(expected_n) +
                                " rows, got " + std::to_string(arr.shape(0)));
  }
  return arr.data();
}

inline AtomMatrix atom_matrix_from_numpy(const NpF64 &arr) {
  if (arr.ndim() != 2 || arr.shape(1) != 3) {
    throw std::invalid_argument(
        "expected float64 array of shape (n, 3), got ndim=" +
        std::to_string(arr.ndim()));
  }
  const auto n = static_cast<Eigen::Index>(arr.shape(0));
  AtomMatrix m(n, 3);
  const double *p = arr.data();
  std::memcpy(m.data(), p, static_cast<size_t>(n) * 3u * sizeof(double));
  return m;
}

inline Matrix3d matrix3_from_numpy(const NpF64 &arr) {
  if (arr.ndim() != 2 || arr.shape(0) != 3 || arr.shape(1) != 3) {
    throw std::invalid_argument("expected float64 array of shape (3, 3)");
  }
  Matrix3d m;
  std::copy(arr.data(), arr.data() + 9, m.data());
  return m;
}

inline VectorXd vector_from_numpy(const NpF64 &arr) {
  if (arr.ndim() != 1) {
    throw std::invalid_argument("expected 1-d float64 array");
  }
  const auto n = static_cast<Eigen::Index>(arr.shape(0));
  VectorXd v(n);
  std::copy(arr.data(), arr.data() + n, v.data());
  return v;
}

inline VectorXi vectori_from_numpy_i64(const NpI64 &arr) {
  if (arr.ndim() != 1) {
    throw std::invalid_argument("expected 1-d int64 array");
  }
  const auto n = static_cast<Eigen::Index>(arr.shape(0));
  VectorXi v(n);
  for (Eigen::Index i = 0; i < n; ++i) {
    v(i) = static_cast<int>(arr.data()[i]);
  }
  return v;
}

// Own a heap buffer so the ndarray remains valid after return.
// Call sites must use nb::rv_policy::move (or automatic): property getters
// default to reference_internal, which conflicts with a capsule owner.
template <typename EigenT>
inline nb::ndarray<nb::numpy, double, nb::c_contig>
matrix_to_numpy(const EigenT &m) {
  const size_t rows = static_cast<size_t>(m.rows());
  const size_t cols = static_cast<size_t>(m.cols());
  double *buf = new double[rows * cols];
  std::copy(m.data(), m.data() + rows * cols, buf);
  // capsule deletes the buffer
  nb::capsule owner(
      buf, [](void *p) noexcept { delete[] static_cast<double *>(p); });
  return nb::ndarray<nb::numpy, double, nb::c_contig>(buf, {rows, cols}, owner);
}

template <typename Derived>
inline nb::ndarray<nb::numpy, double, nb::c_contig>
vector_to_numpy(const Eigen::MatrixBase<Derived> &v) {
  const size_t n = static_cast<size_t>(v.size());
  double *buf = new double[n];
  for (size_t i = 0; i < n; ++i) {
    buf[i] = static_cast<double>(v(static_cast<Eigen::Index>(i)));
  }
  nb::capsule owner(
      buf, [](void *p) noexcept { delete[] static_cast<double *>(p); });
  return nb::ndarray<nb::numpy, double, nb::c_contig>(buf, {n}, owner);
}

/// int64 is the canonical width for _core integer arrays: it round trips
/// through vectori_from_numpy_i64 and matches the numpy default integer.
/// VectorXi storage is int32, so the two differ by a narrowing copy.
inline nb::ndarray<nb::numpy, int64_t, nb::c_contig>
vectori_to_numpy(const VectorXi &v) {
  const size_t n = static_cast<size_t>(v.size());
  int64_t *buf = new int64_t[n];
  for (size_t i = 0; i < n; ++i) {
    buf[i] = static_cast<int64_t>(v(static_cast<Eigen::Index>(i)));
  }
  nb::capsule owner(
      buf, [](void *p) noexcept { delete[] static_cast<int64_t *>(p); });
  return nb::ndarray<nb::numpy, int64_t, nb::c_contig>(buf, {n}, owner);
}

} // namespace eonc::pybind
