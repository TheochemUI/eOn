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
#include "LAMMPSPot.h"
#include "LammpsBundle.h"
#include "LammpsLoader.h"

#include <cstring>
#include <filesystem>
#include <format>
#include <fstream>
#include <map>
#include <string>

#ifdef EONMPI
#define LAMMPS_LIB_MPI
#endif

LAMMPSPot::LAMMPSPot(const Parameters &p)
    : Potential(p),
      lammpsThr{p.potential_options.LAMMPSThreads}
#ifdef EONMPI
      ,
      mpiComm{p.potential_options.MPIClientComm}
#endif
{
  // Fail fast if LAMMPS library not available.
  eonc::LammpsLoader::instance().require_loaded();

  // Two operating modes:
  //   bundle mode  -- lammps_options.bundle_path is set. LAMMPSBundle
  //                   extracts the .eonlpb tarball-equivalent into a
  //                   private scratch dir and we point liblammps at
  //                   that dir via "shell cd <scratch>" so every
  //                   pair_coeff / include / read_data / shared
  //                   plugin .so resolves there. The eonclient CWD
  //                   becomes irrelevant for LAMMPS.
  //   legacy mode  -- bundle_path empty. We require in.lammps in CWD
  //                   so liblammps's relative-path lookups resolve
  //                   against eonclient's CWD; pair-coeff files and
  //                   any other LAMMPS-side file references must
  //                   live there too.
  const auto &bundle_path = p.potential_options.LAMMPSBundlePath;
  if (!bundle_path.empty()) {
    auto bundle = eonc::LAMMPSBundle::open(bundle_path);
    m_lammps_workdir = bundle.extract();
    m_owns_workdir = true;
  } else {
    m_lammps_workdir = std::filesystem::current_path();
    m_owns_workdir = false;
    if (!std::filesystem::exists(m_lammps_workdir / "in.lammps")) {
      auto cwd = m_lammps_workdir.string();
      EONC_LOG_ERROR("[LAMMPS] in.lammps not found in {}", cwd);
      eonc::log::get()->flush_log();
      throw std::runtime_error(
          "LAMMPSPot: in.lammps not found in eonclient CWD (" + cwd +
          "). Either (a) put in.lammps and every file it references "
          "(pair_coeff data, custom pair_style .so plugins, KIM tables, "
          "etc.) next to the eonclient process, or (b) pack them into a "
          "single .eonlpb bundle and pass it via [Potential] "
          "lammps_bundle = path/to/bundle.eonlpb -- liblammps then reads "
          "everything from a private scratch dir, no CWD coupling.");
    }
  }

  // Detect units from in.lammps: look for "#!units real" marker.
  realunits = false;
  std::ifstream infile(m_lammps_workdir / "in.lammps");
  std::string line;
  while (std::getline(infile, line)) {
    if (line == "#!units real") {
      realunits = true;
      break;
    }
  }
}

LAMMPSPot::~LAMMPSPot() {
  cleanMemory();
  if (m_owns_workdir && !m_lammps_workdir.empty()) {
    std::error_code ec;
    std::filesystem::remove_all(m_lammps_workdir, ec);
  }
}

void LAMMPSPot::cleanMemory() {
  if (LAMMPSObj != nullptr) {
    eonc::LammpsLoader::instance().close(LAMMPSObj);
    LAMMPSObj = nullptr;
  }
}

void LAMMPSPot::force(long N, const double *R, const int *atomicNrs, double *F,
                      double *U, double *variance, const double *box) {
  variance = nullptr;
  auto &lmp = eonc::LammpsLoader::instance();

  bool newLammps = false;
  for (int i = 0; i < 9; i++) {
    if (oldBox[i] != box[i])
      newLammps = true;
  }
  if (numberOfAtoms != N)
    newLammps = true;
  if (newLammps) {
    makeNewLAMMPS(N, R, atomicNrs, box);
  }
  if (!LAMMPSObj) {
    throw std::runtime_error("Should have a LAMMPS instance by now");
  }

  lmp.scatter_atoms(LAMMPSObj, "x", 1, 3, const_cast<double *>(R));
  lmp.command(LAMMPSObj, "run 1 pre no post no");

  auto *pe =
      static_cast<double *>(lmp.extract_variable(LAMMPSObj, "pe", nullptr));
  *U = *pe;
  free(pe);

  auto *fx =
      static_cast<double *>(lmp.extract_variable(LAMMPSObj, "fx", "all"));
  auto *fy =
      static_cast<double *>(lmp.extract_variable(LAMMPSObj, "fy", "all"));
  auto *fz =
      static_cast<double *>(lmp.extract_variable(LAMMPSObj, "fz", "all"));

  for (long i = 0; i < N; i++) {
    F[3 * i + 0] = fx[i];
    F[3 * i + 1] = fy[i];
    F[3 * i + 2] = fz[i];
  }

  // Convert kCal/mol -> eV if LAMMPS is using real units
  if (realunits) {
    constexpr double kcalPerEv = 23.0609;
    *U /= kcalPerEv;
    for (long i = 0; i < 3 * N; i++) {
      F[i] /= kcalPerEv;
    }
  }

  free(fx);
  free(fy);
  free(fz);
}

void LAMMPSPot::makeNewLAMMPS(long N, const double *R, const int *atomicNrs,
                              const double *box) {
  auto &lmp = eonc::LammpsLoader::instance();

  numberOfAtoms = N;
  std::memcpy(oldBox, box, 9 * sizeof(double));

  if (LAMMPSObj != nullptr) {
    cleanMemory();
  }

  // Map atomic numbers to LAMMPS type indices (1-based)
  std::map<int, int> type_map;
  int ntypes = 0;
  for (long i = 0; i < N; i++) {
    if (type_map.count(atomicNrs[i]) == 0) {
      type_map.insert({atomicNrs[i], ++ntypes});
    }
  }

#ifdef EONMPI
  const char *lmpargv[] = {"liblammps", "-log", "none",    "-echo", "log",
                           "-screen",   "none", "-suffix", "omp"};
  int lmpargc = sizeof(lmpargv) / sizeof(const char *);
  if (!lmp.open_mpi) {
    throw std::runtime_error(
        "LAMMPS library found but lacks MPI support (lammps_open not found).\n"
        "Install an MPI-enabled LAMMPS build.");
  }
  LAMMPSObj =
      lmp.open_mpi(lmpargc, const_cast<char **>(lmpargv), mpiComm, nullptr);
#else
  const char *lmpargv[] = {"liblammps", "-log",    "none", "-echo",
                           "log",       "-screen", "none"};
  int lmpargc = sizeof(lmpargv) / sizeof(const char *);
  LAMMPSObj = lmp.open_no_mpi(lmpargc, const_cast<char **>(lmpargv), nullptr);
#endif

  // Pin liblammps's working directory to m_lammps_workdir so every
  // pair_coeff / include / read_data / shared-plugin .so reference
  // inside in.lammps resolves there. In bundle mode this is the
  // extracted scratch dir; in legacy mode it's eonclient's CWD (a
  // no-op LAMMPS-side, but harmless).
  {
    std::string shell_cd =
        std::format("shell cd {}", m_lammps_workdir.string());
    lmp.command(LAMMPSObj, shell_cd.c_str());
  }

  if (lammpsThr > 0) {
    std::string cmd = std::format("package omp {} force/neigh", lammpsThr);
    lmp.command(LAMMPSObj, cmd.c_str());
  }

  if (realunits) {
    lmp.command(LAMMPSObj, "units real");
  } else {
    lmp.command(LAMMPSObj, "units metal");
  }

  lmp.command(LAMMPSObj, "atom_style charge");
  lmp.command(LAMMPSObj, "atom_modify map array sort 0 0");
  lmp.command(LAMMPSObj, "neigh_modify delay 1");

  // Define periodic cell (prism for non-orthorhombic)
  std::string region_cmd =
      std::format("region cell prism 0 {} 0 {} 0 {} {} {} {} units box", box[0],
                  box[4], box[8], box[3], box[6], box[7]);
  lmp.command(LAMMPSObj, region_cmd.c_str());

  std::string create_box_cmd = std::format("create_box {} cell", ntypes);
  lmp.command(LAMMPSObj, create_box_cmd.c_str());

  // Initialize atoms
  for (long i = 0; i < N; i++) {
    std::string atom_cmd =
        std::format("create_atoms {} single {} {} {} units box",
                    type_map[atomicNrs[i]], 0.0, 0.0, 0.0);
    lmp.command(LAMMPSObj, atom_cmd.c_str());
  }

  lmp.command(LAMMPSObj, "mass * 1.0");

  // Read in.lammps from the pinned workdir. The "shell cd" above made
  // liblammps's CWD = m_lammps_workdir, so a relative "in.lammps"
  // here resolves against the bundle scratch dir (or eonclient CWD
  // in legacy mode).
  lmp.file(LAMMPSObj, "in.lammps");

  // lammps_file logs syntax / pair_coeff / pair_style errors to its
  // own log and returns void. lammps_has_error (LAMMPS >= 3Mar2020,
  // null on older builds) is the only way to surface them. Without
  // this check the next lammps_command("run 1 ...") would dereference
  // a null pair_style and segfault.
  if (lmp.has_error && lmp.has_error(LAMMPSObj)) {
    char buf[1024]{};
    if (lmp.get_last_error_message) {
      lmp.get_last_error_message(LAMMPSObj, buf, sizeof(buf));
    }
    EONC_LOG_ERROR("[LAMMPS] error after reading in.lammps from {}: {}",
                   m_lammps_workdir.string(), buf);
    eonc::log::get()->flush_log();
    throw std::runtime_error(
        std::string("LAMMPSPot: liblammps reported an error after sourcing ") +
        "in.lammps (workdir=" + m_lammps_workdir.string() +
        "). LAMMPS message: " + buf);
  }

  // Define variables for force/energy extraction
  lmp.command(LAMMPSObj, "variable fx atom fx");
  lmp.command(LAMMPSObj, "variable fy atom fy");
  lmp.command(LAMMPSObj, "variable fz atom fz");
  lmp.command(LAMMPSObj, "variable pe equal pe");
}
