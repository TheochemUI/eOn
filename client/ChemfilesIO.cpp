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
#include "eon/ChemfilesIO.h"
#include "eon/EonLogger.h"
#include "eon/Matter.h"

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <format>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <string>

#ifdef WITH_CHEMFILES
#include <chemfiles.hpp>
#endif

namespace eonc::io {
namespace {

namespace fs = std::filesystem;

std::string lower(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  return s;
}

std::mutex &format_mutex() {
  static std::mutex mu;
  return mu;
}
std::string &input_fmt() {
  static std::string s;
  return s;
}
thread_local std::string g_last_cell_warning;

struct ConfigCell {
  std::vector<double> values; // empty, 3 or 9
  bool pbc{true};
  double vacuum{10.0};
};
ConfigCell &config_cell() {
  static ConfigCell c;
  return c;
}

/// Cell from [Main] cell, rows are lattice vectors.
Matrix3d matrix_from_values(const std::vector<double> &v) {
  Matrix3d box = Matrix3d::Zero();
  if (v.size() == 3) {
    box(0, 0) = v[0];
    box(1, 1) = v[1];
    box(2, 2) = v[2];
  } else {
    for (int i = 0; i < 9; ++i) {
      box(i / 3, i % 3) = v[static_cast<size_t>(i)];
    }
  }
  return box;
}

#ifdef WITH_CHEMFILES
/// Map an eOn format word to a chemfiles format name ("" = guess).
std::string chemfiles_format(const std::string &fmt) {
  const auto f = lower(fmt);
  if (f.empty()) {
    return "";
  }
  if (f == "xyz" || f == "extxyz") {
    return "XYZ";
  }
  if (f == "pdb") {
    return "PDB";
  }
  if (f == "cif") {
    return "CIF";
  }
  if (f == "gro") {
    return "GRO";
  }
  if (f == "data" || f == "lammps" || f == "lammps-data" || f == "lmp") {
    return "LAMMPS Data";
  }
  if (f == "mol2") {
    return "MOL2";
  }
  if (f == "sdf") {
    return "SDF";
  }
  if (f == "poscar" || f == "vasp") {
    return "POSCAR";
  }
  return fmt; // a chemfiles format name, passed through
}

/// Extensions chemfiles cannot guess but eOn users write.
std::string format_for_path(const std::string &filename,
                            const std::string &format) {
  if (!format.empty()) {
    return chemfiles_format(format);
  }
  const auto ext = lower(fs::path(filename).extension().string());
  if (ext == ".extxyz") {
    return "XYZ";
  }
  if (ext == ".data" || ext == ".lmp") {
    return "LAMMPS Data";
  }
  const auto name = lower(fs::path(filename).filename().string());
  if (name == "poscar" || name == "contcar" || ext == ".vasp") {
    return "POSCAR";
  }
  return "";
}

/// Read a fixed mask from one property. `frozen_when_true` is false for
/// move_mask, where true means the atom moves.
bool mask_from_property(const chemfiles::Property &p, bool frozen_when_true,
                        std::array<bool, 3> &mask) {
  auto flip = [&](bool v) { return frozen_when_true ? v : !v; };
  switch (p.kind()) {
  case chemfiles::Property::BOOL: {
    const bool v = flip(p.as_bool());
    mask = {v, v, v};
    return true;
  }
  case chemfiles::Property::DOUBLE: {
    const bool v = flip(p.as_double() != 0.0);
    mask = {v, v, v};
    return true;
  }
  case chemfiles::Property::VECTOR3D: {
    const auto v = p.as_vector3d();
    mask = {flip(v[0] != 0.0), flip(v[1] != 0.0), flip(v[2] != 0.0)};
    return true;
  }
  default:
    return false;
  }
}

void quiet_chemfiles() {
  static std::once_flag once;
  std::call_once(once, [] {
    chemfiles::set_warning_callback([](const std::string &msg) {
      EONC_LOG_DEBUG("chemfiles: {}", msg);
    });
  });
}
#endif

} // namespace

bool chemfiles_available() noexcept {
#ifdef WITH_CHEMFILES
  return true;
#else
  return false;
#endif
}

bool is_readcon_path(std::string_view filename) {
  const auto name = lower(fs::path(std::string(filename)).filename().string());
  for (const char *ext : {".con", ".convel"}) {
    const std::string e(ext);
    if (name.ends_with(e) || name.ends_with(e + ".gz") ||
        name.ends_with(e + ".zst")) {
      return true;
    }
  }
  return false;
}

bool is_chemfiles_path(std::string_view filename) {
  const fs::path path{std::string(filename)};
  const auto name = lower(path.filename().string());
  if (name == "poscar" || name == "contcar") {
    return true;
  }
  static const char *const kExt[] = {
      ".xyz", ".extxyz", ".pdb", ".cif",  ".mmcif", ".gro", ".data",
      ".lmp", ".mol2",   ".sdf", ".vasp", ".xtc",   ".trr", ".dcd",
      ".nc",  ".arc",    ".cml", ".cssr", ".mmtf",  ".tng"};
  const auto ext = lower(path.extension().string());
  return std::find(std::begin(kExt), std::end(kExt), ext) != std::end(kExt);
}

void set_structure_input_format(std::string fmt) {
  fmt = lower(std::move(fmt));
  const std::lock_guard<std::mutex> g(format_mutex());
  input_fmt() = (fmt == "con") ? std::string() : fmt;
}

const std::string &structure_input_format() noexcept { return input_fmt(); }

std::vector<double> parse_cell_values(const std::string &text) {
  std::vector<double> out;
  std::istringstream in(text);
  std::string tok;
  while (in >> tok) {
    try {
      size_t used = 0;
      out.push_back(std::stod(tok, &used));
      if (used != tok.size()) {
        throw std::invalid_argument(tok);
      }
    } catch (const std::exception &) {
      throw std::invalid_argument("[Main] cell: '" + tok +
                                  "' is not a number");
    }
  }
  if (!out.empty() && out.size() != 3 && out.size() != 9) {
    throw std::invalid_argument(
        "[Main] cell takes 3 numbers (a b c) or 9 (a 3x3 lattice), got " +
        std::to_string(out.size()));
  }
  return out;
}

void set_config_cell(const std::vector<double> &values, bool pbc) {
  if (!values.empty() && values.size() != 3 && values.size() != 9) {
    throw std::invalid_argument("[Main] cell takes 3 or 9 numbers");
  }
  const std::lock_guard<std::mutex> g(format_mutex());
  config_cell().values = values;
  config_cell().pbc = pbc;
}

void set_default_vacuum(double angstrom) {
  if (!(angstrom >= 0.0)) {
    throw std::invalid_argument("[Main] default_vacuum must be >= 0");
  }
  const std::lock_guard<std::mutex> g(format_mutex());
  config_cell().vacuum = angstrom;
}

double default_vacuum() noexcept { return config_cell().vacuum; }

const std::string &last_cell_warning() noexcept { return g_last_cell_warning; }

IoStatus inp2matter(Matter &m, const std::string &filename) {
  if (is_chemfiles_path(filename)) {
    return chemfiles2matter(m, filename);
  }
  const auto name = lower(fs::path(filename).filename().string());
  if (name.ends_with(".convel") || name.ends_with(".convel.gz") ||
      name.ends_with(".convel.zst")) {
    return convel2matter(m, filename);
  }
  return con2matter(m, filename);
}

#ifdef WITH_CHEMFILES

IoStatus chemfiles2matter(Matter &m, const std::string &filename,
                          const std::string &format) {
  quiet_chemfiles();
  try {
    chemfiles::Trajectory traj(filename, 'r',
                               format_for_path(filename, format));
    if (traj.nsteps() == 0) {
      EONC_LOG_ERROR("{} holds no frames", filename);
      return IoStatus::ReadError;
    }
    const chemfiles::Frame frame = traj.read();

    g_last_cell_warning.clear();
    const auto &cell = frame.cell();
    const bool file_has_cell = cell.shape() != chemfiles::UnitCell::INFINITE;
    const auto n = static_cast<Eigen::Index>(frame.size());
    m.resize(static_cast<long>(n));

    AtomMatrix positions(n, 3);
    VectorXd masses(n);
    VectorXi atomic_nrs(n);
    const auto &pos = frame.positions();
    for (Eigen::Index i = 0; i < n; ++i) {
      const auto &atom = frame[static_cast<size_t>(i)];
      const auto z = atom.atomic_number();
      if (!z) {
        EONC_LOG_ERROR("{}: atom {} has type '{}', which is not an element",
                       filename, i, atom.type());
        return IoStatus::InvalidArgument;
      }
      atomic_nrs(i) = static_cast<int>(*z);
      masses(i) = atom.mass();
      for (int k = 0; k < 3; ++k) {
        positions(i, k) = pos[static_cast<size_t>(i)][k];
      }
      std::array<bool, 3> mask{false, false, false};
      if (auto p = atom.get("fixed")) {
        mask_from_property(*p, /*frozen_when_true=*/true, mask);
      } else if (auto p = atom.get("move_mask")) {
        mask_from_property(*p, /*frozen_when_true=*/false, mask);
      }
      m.setFixedMask(static_cast<long>(i), mask);
      m.setAtomIndex(static_cast<long>(i), static_cast<std::int64_t>(i));
    }
    // chemfiles keeps file order, so the file-to-Matter map is the identity.
    std::vector<long> identity(static_cast<size_t>(n));
    for (Eigen::Index i = 0; i < n; ++i) {
      identity[static_cast<size_t>(i)] = static_cast<long>(i);
    }
    m.setFileToMatter(std::move(identity));

    Matrix3d box = Matrix3d::Zero();
    bool periodic = true;
    const ConfigCell cfg = [] {
      const std::lock_guard<std::mutex> g(format_mutex());
      return config_cell();
    }();
    if (file_has_cell) {
      // chemfiles stores cell vectors as columns; Matter stores rows.
      const auto cm = cell.matrix();
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          box(i, j) = cm[j][i];
        }
      }
    } else if (!cfg.values.empty()) {
      box = matrix_from_values(cfg.values);
      periodic = cfg.pbc;
    } else {
      // Bounding box plus vacuum, shifted so every atom sits inside.
      Eigen::RowVector3d lo = Eigen::RowVector3d::Zero();
      Eigen::RowVector3d hi = Eigen::RowVector3d::Zero();
      if (n > 0) {
        lo = positions.colwise().minCoeff();
        hi = positions.colwise().maxCoeff();
      }
      for (int k = 0; k < 3; ++k) {
        box(k, k) = (hi(k) - lo(k)) + 2.0 * cfg.vacuum;
        // A flat or single-atom axis still gets a positive length.
        if (box(k, k) <= 0.0) {
          box(k, k) = 1.0;
        }
      }
      positions.rowwise() -= (lo - Eigen::RowVector3d::Constant(cfg.vacuum));
      periodic = false;
      g_last_cell_warning = std::format(
          "WARNING: {} has no unit cell and [Main] cell is not set. Using a "
          "non-periodic orthorhombic box of {:.3f} x {:.3f} x {:.3f} A "
          "(atom bounding box plus {} A vacuum per side, [Main] "
          "default_vacuum); atoms were shifted into it. Set [Main] cell and "
          "[Main] pbc, or add a cell to the file, for a periodic system.",
          filename, box(0, 0), box(1, 1), box(2, 2), cfg.vacuum);
      EONC_LOG_WARNING("{}", g_last_cell_warning);
    }
    m.setCell(box);
    if (!periodic) {
      m.setPeriodic(false); // never force PBC onto an isolated-layout pot
    }
    m.setMasses(masses);
    m.setAtomicNrs(atomic_nrs);
    m.setPositions(positions);

    // Frame::velocities() is set by formats that store velocities as a
    // trajectory array. The XYZ reader does not (FormatMetadata::velocities
    // is false): an extxyz `velo` column is an atom property of kind
    // VECTOR3D, so read that when the frame has no velocity array.
    AtomMatrix v(n, 3);
    bool have_velocities = false;
    if (const auto vel = frame.velocities()) {
      have_velocities = true;
      for (Eigen::Index i = 0; i < n; ++i) {
        for (int k = 0; k < 3; ++k) {
          v(i, k) = (*vel)[static_cast<size_t>(i)][k];
        }
      }
    } else {
      v.setZero();
      for (Eigen::Index i = 0; i < n; ++i) {
        const auto &atom = frame[static_cast<size_t>(i)];
        const auto prop = atom.get("velo");
        if (!prop || prop->kind() != chemfiles::Property::VECTOR3D) {
          continue;
        }
        have_velocities = true;
        const auto vec = prop->as_vector3d();
        for (int k = 0; k < 3; ++k) {
          v(i, k) = vec[static_cast<size_t>(k)];
        }
      }
    }
    if (have_velocities) {
      m.setVelocities(v);
    }
    // resize/setPositions leave the potential dirty; no file forces here.
    return IoStatus::Ok;
  } catch (const std::exception &e) {
    EONC_LOG_ERROR("Failed to read {} with chemfiles: {}", filename, e.what());
    return IoStatus::ReadError;
  }
}

#else // !WITH_CHEMFILES

IoStatus chemfiles2matter(Matter &, const std::string &filename,
                          const std::string &) {
  EONC_LOG_ERROR("Cannot read {}: this eonclient was built without chemfiles. "
                 "Install chemfiles-lib and reconfigure with "
                 "-Dwith_chemfiles=enabled, or convert the file to .con.",
                 filename);
  return IoStatus::InvalidArgument;
}

#endif

} // namespace eonc::io
