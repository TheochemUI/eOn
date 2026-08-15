#pragma once
/** Shared Matter helpers for pyeonclient bindings only (not on Matter itself).
 *
 * Contract: algorithms never mutate the caller's Matter unless
 * inplace=True. Default is always non-mutating (copy-then-run).
 *
 * Buffer access uses public Matter getters/setters; ASE/Python cache state
 * lives in the binding layer (AseCalcPotential), not on Matter.
 */
#include "eon/ConFileIO.h"
#include "eon/Eigen.h"
#include "eon/Matter.h"

#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>

#include <filesystem>
#include <format>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace eonc::pybind {
namespace nb = nanobind;

/// Working Matter: same object if inplace, otherwise deep copy.
inline std::shared_ptr<eonc::Matter>
matter_work(std::shared_ptr<eonc::Matter> matter, bool inplace) {
  if (inplace)
    return matter;
  return std::make_shared<eonc::Matter>(*matter);
}

/// Matter stores a non-owning `const Parameters *` (Matter.h) and every copy
/// shares it, so a derived Matter handed to Python must hold a reference to
/// something that owns those Parameters. `nb::keep_alive<0, N>` covers a
/// return value that *is* the object; it cannot reach inside a returned tuple
/// or list, so those sites tie their elements through here.
inline void tie_lifetime(nb::handle nurse, nb::handle patient) {
  if (!nurse.is_valid() || !patient.is_valid())
    return;
  if (nurse.is_none() || patient.is_none() || nurse.ptr() == patient.ptr())
    return;
  nb::detail::keep_alive(nurse.ptr(), patient.ptr());
}

/// Python object already wrapping `m`, or an invalid handle when Python has
/// never seen it.
inline nb::object matter_object(const eonc::Matter &m) { return nb::find(m); }

/// Cast a freshly built path to a list, tying each frame to the object whose
/// Parameters the frames point at.
inline nb::list matter_path_to_python(std::vector<eonc::Matter> path,
                                      nb::handle owner) {
  nb::list out;
  for (auto &frame : path) {
    nb::object obj = nb::cast(std::move(frame));
    tie_lifetime(obj, owner);
    out.append(obj);
  }
  return out;
}

/// Write C++ ConFrames to a short-lived temp .con and load as Python
/// readcon.ConFrame list (same stamps as disk). Caller never sees the path.
/// Takes const ref: readcon::ConFrame is move-only (copy deleted).
inline nb::list
con_frames_to_python(const std::vector<readcon::ConFrame> &frames) {
  if (frames.empty()) {
    return nb::list();
  }
  namespace fs = std::filesystem;
  const auto tmp =
      fs::temp_directory_path() /
      std::format("eon_path_frames_{}_{}.con",
                  static_cast<unsigned>(std::random_device{}() % 1000000u),
                  reinterpret_cast<uintptr_t>(frames.data()));
  if (!eonc::io::io_ok(eonc::io::writeConFrames(tmp.string(), frames))) {
    throw std::runtime_error("path_frames: failed to serialize ConFrames");
  }
  try {
    nb::object readcon = nb::module_::import_("readcon");
    nb::object py_frames = readcon.attr("read_con")(tmp.string());
    fs::remove(tmp);
    return nb::cast<nb::list>(py_frames);
  } catch (...) {
    std::error_code ec;
    fs::remove(tmp, ec);
    throw;
  }
}

/// Fold algorithm result into the caller's Matter when inplace=True.
inline std::shared_ptr<eonc::Matter>
matter_result(std::shared_ptr<eonc::Matter> original,
              std::shared_ptr<eonc::Matter> result, bool inplace) {
  if (!result)
    return result;
  if (!inplace)
    return result;
  if (result.get() != original.get())
    *original = *result;
  return original;
}

// Binding-local zero-copy views / bulk sets via public Matter API.

inline double *matter_positions_ptr(eonc::Matter &m) {
  return const_cast<double *>(m.getPositions().data());
}
inline const double *matter_positions_ptr(const eonc::Matter &m) {
  return m.getPositions().data();
}
inline double *matter_forces_ptr(eonc::Matter &m) {
  return const_cast<double *>(m.getForces().data());
}
inline double *matter_forces_raw_ptr(eonc::Matter &m) {
  return const_cast<double *>(m.getForcesRaw().data());
}

inline void matter_set_positions_buf(eonc::Matter &m, const double *xyz,
                                     long n) {
  if (n != m.numberOfAtoms()) {
    throw std::invalid_argument("positions n != n_atoms");
  }
  m.setPositions(AtomMatrix::Map(const_cast<double *>(xyz), n, 3));
}

inline void matter_set_cell_buf(eonc::Matter &m, const double *c33) {
  m.setCell(Matrix3d::Map(const_cast<double *>(c33)));
}

inline void matter_set_velocities_buf(eonc::Matter &m, const double *v,
                                      long n) {
  if (n != m.numberOfAtoms()) {
    throw std::invalid_argument("velocities n != n_atoms");
  }
  m.setVelocities(AtomMatrix::Map(const_cast<double *>(v), n, 3));
}

inline void matter_set_masses_buf(eonc::Matter &m, const double *mass, long n) {
  if (n != m.numberOfAtoms()) {
    throw std::invalid_argument("masses n != n_atoms");
  }
  m.setMasses(
      Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, 1>>(mass, n));
}

inline void matter_set_z_buf(eonc::Matter &m, const int *z, long n) {
  if (n != m.numberOfAtoms()) {
    throw std::invalid_argument("atomic_numbers n != n_atoms");
  }
  m.setAtomicNrs(VectorXi::Map(const_cast<int *>(z), n));
}

inline void matter_set_fixed_buf(eonc::Matter &m, const int *fx, long n) {
  if (n != m.numberOfAtoms()) {
    throw std::invalid_argument("fixed n != n_atoms");
  }
  for (long i = 0; i < n; ++i) {
    m.setFixed(i, fx[i] ? 1 : 0);
  }
}

} // namespace eonc::pybind
