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
#include "eon/ParametersINI.h"
#include "eon/BaseStructures.h"
#include "eon/ConFileIO.h"
#include "eon/EpiCenters.h"
#include "eon/HelperFunctions.h"
#include "eon/Parameters.h"
#include "magic_enum/magic_enum.hpp"

#include <INIReader.h>

#include <cerrno>
#include <cmath>
#include <ctime>
#include <sstream>
#include <stdexcept>
#include <string>

#include "eon/EonLogger.h"

namespace {
std::string toLowerCase(std::string s) {
  for (std::string::size_type i = 0; i < s.length(); ++i) {
    s[i] = tolower(s[i]);
  }
  return s;
}
} // namespace

namespace eonc::config {

int load_ini(INIReader &ini, Parameters &params) {
  int error = 0;

  // [Main] //

  ParametersLoadAccess::main_options(params).job =
      magic_enum::enum_cast<JobType>(ini.Get("Main", "job", ""),
                                     magic_enum::case_insensitive)
          .value_or(JobType::Unknown);
  ParametersLoadAccess::main_options(params).temperature =
      ini.GetReal("Main", "temperature",
                  ParametersLoadAccess::main_options(params).temperature);
  ParametersLoadAccess::main_options(params).randomSeed =
      ini.GetInteger("Main", "random_seed",
                     ParametersLoadAccess::main_options(params).randomSeed);
  ParametersLoadAccess::main_options(params).checkpoint =
      ini.GetBoolean("Main", "checkpoint",
                     ParametersLoadAccess::main_options(params).checkpoint);
  ParametersLoadAccess::main_options(params).quiet = ini.GetBoolean(
      "Main", "quiet", ParametersLoadAccess::main_options(params).quiet);
  ParametersLoadAccess::main_options(params).writeLog = ini.GetBoolean(
      "Main", "write_log", ParametersLoadAccess::main_options(params).writeLog);
  ParametersLoadAccess::main_options(params).writeConForces =
      ini.GetBoolean("Main", "write_con_forces",
                     ParametersLoadAccess::main_options(params).writeConForces);
  // Process-wide io mode; same pattern as the RNG seeding below.
  eonc::io::set_write_con_forces(
      ParametersLoadAccess::main_options(params).writeConForces);
  ParametersLoadAccess::main_options(params).finiteDifference =
      ini.GetReal("Main", "finite_difference",
                  ParametersLoadAccess::main_options(params).finiteDifference);
  // Initialize random generator
  if (ParametersLoadAccess::main_options(params).randomSeed < 0) {
    unsigned i = static_cast<unsigned>(std::time(nullptr));
    ParametersLoadAccess::main_options(params).randomSeed = i;
    eonc::rng::random(i);
  } else {
    eonc::rng::random(ParametersLoadAccess::main_options(params).randomSeed);
  }
  ParametersLoadAccess::main_options(params).maxForceCalls =
      ini.GetInteger("Main", "max_force_calls",
                     ParametersLoadAccess::main_options(params).maxForceCalls);
  ParametersLoadAccess::main_options(params).removeNetForce =
      ini.GetBoolean("Main", "remove_net_force",
                     ParametersLoadAccess::main_options(params).removeNetForce);
  ParametersLoadAccess::main_options(params).parallel = ini.GetBoolean(
      "Main", "parallel", ParametersLoadAccess::main_options(params).parallel);

  // [Potential] //

  std::string potTok = ini.Get("Potential", "potential", "");
  // Schema / old configs: ase_nwcem is ASE_NWCHEM; socket_nwchem is
  // SocketNWChem (magic_enum already matches the latter).
  if (potTok == "ase_nwcem" || potTok == "ASE_NWCEM") {
    potTok = "ase_nwchem";
  }
  ParametersLoadAccess::potential_options(params).potential =
      magic_enum::enum_cast<PotType>(potTok, magic_enum::case_insensitive)
          .value_or(PotType::UNKNOWN);
  ParametersLoadAccess::potential_options(params).MPIPollPeriod = ini.GetReal(
      "Potential", "mpi_poll_period",
      ParametersLoadAccess::potential_options(params).MPIPollPeriod);
  ParametersLoadAccess::potential_options(params).LAMMPSLogging =
      ini.GetBoolean(
          "Potential", "lammps_logging",
          ParametersLoadAccess::potential_options(params).LAMMPSLogging);
  ParametersLoadAccess::potential_options(params).LAMMPSThreads =
      static_cast<int>(ini.GetInteger(
          "Potential", "lammps_threads",
          ParametersLoadAccess::potential_options(params).LAMMPSThreads));
  ParametersLoadAccess::potential_options(params).EMTRasmussen = ini.GetBoolean(
      "Potential", "emt_rasmussen",
      ParametersLoadAccess::potential_options(params).EMTRasmussen);
  ParametersLoadAccess::potential_options(params).extPotPath =
      ini.Get("Potential", "ext_pot_path",
              ParametersLoadAccess::potential_options(params).extPotPath);
  ParametersLoadAccess::potential_options(params).potentialsPath =
      ini.Get("Potential", "potentials_path",
              ParametersLoadAccess::potential_options(params).potentialsPath);

  if (params.potential_options().potential == PotType::MPI ||
      params.potential_options().potential == PotType::VASP) {
    ParametersLoadAccess::potential_options(params).LogPotential = true;
  } else {
    ParametersLoadAccess::potential_options(params).LogPotential = false;
  }
  ParametersLoadAccess::potential_options(params).LogPotential = ini.GetBoolean(
      "Potential", "log_potential",
      ParametersLoadAccess::potential_options(params).LogPotential);
  ParametersLoadAccess::potential_options(params).thread_safe = ini.GetBoolean(
      "Potential", "thread_safe",
      ParametersLoadAccess::potential_options(params).thread_safe);

  // [AMS]
  if (params.potential_options().potential == PotType::AMS) {
    ParametersLoadAccess::ams_options(params).engine = ini.Get(
        "AMS", "engine", ParametersLoadAccess::ams_options(params).engine);
    ParametersLoadAccess::ams_options(params).forcefield =
        ini.Get("AMS", "forcefield",
                ParametersLoadAccess::ams_options(params).forcefield);
    ParametersLoadAccess::ams_options(params).resources =
        ini.Get("AMS", "resources",
                ParametersLoadAccess::ams_options(params).resources);
    ParametersLoadAccess::ams_options(params).model = ini.Get(
        "AMS", "model", ParametersLoadAccess::ams_options(params).model);
    ParametersLoadAccess::ams_options(params).xc =
        ini.Get("AMS", "xc", ParametersLoadAccess::ams_options(params).xc);
    ParametersLoadAccess::ams_options(params).basis = ini.Get(
        "AMS", "basis", ParametersLoadAccess::ams_options(params).basis);
  }
  // [AMS_IO]
  if (params.potential_options().potential == PotType::AMS_IO) {
    ParametersLoadAccess::ams_options(params).engine = ini.Get(
        "AMS_IO", "engine", ParametersLoadAccess::ams_options(params).engine);
    ParametersLoadAccess::ams_options(params).forcefield =
        ini.Get("AMS_IO", "forcefield",
                ParametersLoadAccess::ams_options(params).forcefield);
    ParametersLoadAccess::ams_options(params).model = ini.Get(
        "AMS_IO", "model", ParametersLoadAccess::ams_options(params).model);
    ParametersLoadAccess::ams_options(params).xc =
        ini.Get("AMS_IO", "xc", ParametersLoadAccess::ams_options(params).xc);
  }
  // [AMS_ENV]
  if (params.potential_options().potential == PotType::AMS_IO ||
      params.potential_options().potential == PotType::AMS) {
    ParametersLoadAccess::ams_options(params).env.amshome =
        ini.Get("AMS_ENV", "amshome",
                ParametersLoadAccess::ams_options(params).env.amshome);
    ParametersLoadAccess::ams_options(params).env.scm_tmpdir =
        ini.Get("AMS_ENV", "scm_tmpdir",
                ParametersLoadAccess::ams_options(params).env.scm_tmpdir);
    ParametersLoadAccess::ams_options(params).env.scmlicense =
        ini.Get("AMS_ENV", "scmlicense",
                ParametersLoadAccess::ams_options(params).env.scmlicense);
    ParametersLoadAccess::ams_options(params).env.scm_pythondir =
        ini.Get("AMS_ENV", "scm_pythondir",
                ParametersLoadAccess::ams_options(params).env.scm_pythondir);
    ParametersLoadAccess::ams_options(params).env.amsbin =
        ini.Get("AMS_ENV", "amsbin",
                ParametersLoadAccess::ams_options(params).env.amsbin);
    ParametersLoadAccess::ams_options(params).env.amsresources =
        ini.Get("AMS_ENV", "amsresources",
                ParametersLoadAccess::ams_options(params).env.amsresources);
  }
  // [XTBPot]
  if (params.potential_options().potential == PotType::XTB) {
    ParametersLoadAccess::xtb_options(params).paramset =
        ini.Get("XTBPot", "paramset",
                ParametersLoadAccess::xtb_options(params).paramset);
    ParametersLoadAccess::xtb_options(params).acc = ini.GetReal(
        "XTBPot", "accuracy", ParametersLoadAccess::xtb_options(params).acc);
    ParametersLoadAccess::xtb_options(params).elec_temperature =
        ini.GetReal("XTBPot", "electronic_temperature",
                    ParametersLoadAccess::xtb_options(params).elec_temperature);
    ParametersLoadAccess::xtb_options(params).maxiter =
        ini.GetInteger("XTBPot", "max_iterations",
                       ParametersLoadAccess::xtb_options(params).maxiter);
    ParametersLoadAccess::xtb_options(params).uhf = ini.GetInteger(
        "XTBPot", "uhf", ParametersLoadAccess::xtb_options(params).uhf);
    ParametersLoadAccess::xtb_options(params).charge = ini.GetReal(
        "XTBPot", "charge", ParametersLoadAccess::xtb_options(params).charge);
  }
  // [ZBLPot]
  if (params.potential_options().potential == PotType::ZBL) {
    ParametersLoadAccess::zbl_options(params).cut_inner =
        ini.GetReal("ZBLPot", "cut_inner",
                    ParametersLoadAccess::zbl_options(params).cut_inner);
    ParametersLoadAccess::zbl_options(params).cut_global =
        ini.GetReal("ZBLPot", "cut_global",
                    ParametersLoadAccess::zbl_options(params).cut_global);
    if (ParametersLoadAccess::zbl_options(params).cut_inner >
        ParametersLoadAccess::zbl_options(params).cut_global) {
      throw std::runtime_error(
          "Switching function must begin before the global cutoff!");
    }
  }
  // [D3Pot] / [D4Pot]: rgpot 3.1 Grimme DFT-D
  if (params.potential_options().potential == PotType::DFTD3 ||
      params.potential_options().potential == PotType::DFTD4) {
    ParametersLoadAccess::dftd_options(params).functional =
        ini.Get("D3Pot", "functional",
                ini.Get("D4Pot", "functional",
                        ParametersLoadAccess::dftd_options(params).functional));
    ParametersLoadAccess::dftd_options(params).atm = ini.GetBoolean(
        "D3Pot", "atm",
        ini.GetBoolean("D4Pot", "atm",
                       ParametersLoadAccess::dftd_options(params).atm));
    ParametersLoadAccess::dftd_options(params).d3_damping =
        ini.Get("D3Pot", "damping",
                ParametersLoadAccess::dftd_options(params).d3_damping);
    ParametersLoadAccess::dftd_options(params).d4_charge =
        ini.GetReal("D4Pot", "charge",
                    ParametersLoadAccess::dftd_options(params).d4_charge);
  }
  if (params.potential_options().potential == PotType::EXPR) {
    ParametersLoadAccess::expr_options(params).expression =
        ini.Get("ExprPot", "expression",
                ParametersLoadAccess::expr_options(params).expression);
    ParametersLoadAccess::expr_options(params).terms = ini.Get(
        "ExprPot", "terms", ParametersLoadAccess::expr_options(params).terms);
  }
  if (params.potential_options().potential == PotType::MOPAC) {
    ParametersLoadAccess::mopac_options(params).charge = static_cast<int>(
        ini.GetInteger("MOPACPot", "charge",
                       ParametersLoadAccess::mopac_options(params).charge));
    ParametersLoadAccess::mopac_options(params).spin = static_cast<int>(
        ini.GetInteger("MOPACPot", "spin",
                       ParametersLoadAccess::mopac_options(params).spin));
    ParametersLoadAccess::mopac_options(params).model = static_cast<int>(
        ini.GetInteger("MOPACPot", "model",
                       ParametersLoadAccess::mopac_options(params).model));
    ParametersLoadAccess::mopac_options(params).engine_path =
        ini.Get("MOPACPot", "engine_path",
                ParametersLoadAccess::mopac_options(params).engine_path);
  }
  // [SocketNWChemPot]
  if (params.potential_options().potential == PotType::SocketNWChem) {
    ParametersLoadAccess::socket_nwchem_options(params).host =
        ini.Get("SocketNWChemPot", "host",
                ParametersLoadAccess::socket_nwchem_options(params).host);
    ParametersLoadAccess::socket_nwchem_options(params).port = ini.GetInteger(
        "SocketNWChemPot", "port",
        ParametersLoadAccess::socket_nwchem_options(params).port);
    ParametersLoadAccess::socket_nwchem_options(params).mem_in_gb =
        ini.GetInteger(
            "SocketNWChemPot", "mem_in_gb",
            ParametersLoadAccess::socket_nwchem_options(params).mem_in_gb);
    ParametersLoadAccess::socket_nwchem_options(params).nwchem_settings =
        ini.Get("SocketNWChemPot", "nwchem_settings",
                ParametersLoadAccess::socket_nwchem_options(params)
                    .nwchem_settings);
    ParametersLoadAccess::socket_nwchem_options(params).unix_socket_path =
        ini.Get("SocketNWChemPot", "unix_socket_path",
                ParametersLoadAccess::socket_nwchem_options(params)
                    .unix_socket_path);
    ParametersLoadAccess::socket_nwchem_options(params).unix_socket_mode =
        ini.GetBoolean("SocketNWChemPot", "unix_socket_mode",
                       ParametersLoadAccess::socket_nwchem_options(params)
                           .unix_socket_mode);
    ParametersLoadAccess::socket_nwchem_options(params).make_template_input =
        ini.GetBoolean("SocketNWChemPot", "make_template_input",
                       ParametersLoadAccess::socket_nwchem_options(params)
                           .make_template_input);
  }

  // [RgpotPot] — in-process NWChemPot/CPMDPot (also accept legacy [RGPot] keys)
  if (params.potential_options().potential == PotType::RGPOT) {
    const char *sec = "RgpotPot";
    // Prefer [RgpotPot]; fall back to [RGPot] field names used by direct-link
    // design
    ParametersLoadAccess::rgpot_options(params).backend = ini.Get(
        sec, "backend", ParametersLoadAccess::rgpot_options(params).backend);
    ParametersLoadAccess::rgpot_options(params).basis =
        ini.Get(sec, "basis",
                ini.Get(sec, "nwchem_basis",
                        ParametersLoadAccess::rgpot_options(params).basis));
    ParametersLoadAccess::rgpot_options(params).theory =
        ini.Get(sec, "theory",
                ini.Get(sec, "nwchem_theory",
                        ParametersLoadAccess::rgpot_options(params).theory));
    ParametersLoadAccess::rgpot_options(params).scf_type =
        ini.Get(sec, "scf_type",
                ini.Get(sec, "nwchem_scf_type",
                        ParametersLoadAccess::rgpot_options(params).scf_type));
    ParametersLoadAccess::rgpot_options(params).functional = ini.Get(
        sec, "functional",
        ini.Get(sec, "cpmd_functional",
                ParametersLoadAccess::rgpot_options(params).functional));
    ParametersLoadAccess::rgpot_options(params).cutoff_ry = ini.GetReal(
        sec, "cutoff_ry",
        ini.GetReal(sec, "cpmd_cut_off_ry",
                    ParametersLoadAccess::rgpot_options(params).cutoff_ry));
    ParametersLoadAccess::rgpot_options(params).charge = ini.GetInteger(
        sec, "charge",
        ini.GetInteger(sec, "nwchem_charge",
                       ParametersLoadAccess::rgpot_options(params).charge));
    ParametersLoadAccess::rgpot_options(params).multiplicity = ini.GetInteger(
        sec, "multiplicity",
        ini.GetInteger(
            sec, "nwchem_multiplicity",
            ParametersLoadAccess::rgpot_options(params).multiplicity));
    ParametersLoadAccess::rgpot_options(params).engine_path =
        ini.Get(sec, "engine_path",
                ParametersLoadAccess::rgpot_options(params).engine_path);
    ParametersLoadAccess::rgpot_options(params).engine_library =
        ini.Get(sec, "engine_library",
                ParametersLoadAccess::rgpot_options(params).engine_library);
    ParametersLoadAccess::rgpot_options(params).engine_root =
        ini.Get(sec, "engine_root",
                ParametersLoadAccess::rgpot_options(params).engine_root);
    ParametersLoadAccess::rgpot_options(params).title = ini.Get(
        sec, "title", ParametersLoadAccess::rgpot_options(params).title);
    ParametersLoadAccess::rgpot_options(params).memory_mb =
        ini.GetInteger(sec, "memory_mb",
                       ParametersLoadAccess::rgpot_options(params).memory_mb);
    ParametersLoadAccess::rgpot_options(params).scratch_dir =
        ini.Get(sec, "scratch_dir",
                ParametersLoadAccess::rgpot_options(params).scratch_dir);
    ParametersLoadAccess::rgpot_options(params).input_block =
        ini.Get(sec, "input_block",
                ParametersLoadAccess::rgpot_options(params).input_block);
    ParametersLoadAccess::rgpot_options(params).model_path =
        ini.Get(sec, "model_path",
                ParametersLoadAccess::rgpot_options(params).model_path);
    ParametersLoadAccess::rgpot_options(params).device = ini.Get(
        sec, "device", ParametersLoadAccess::rgpot_options(params).device);
    ParametersLoadAccess::rgpot_options(params).length_unit =
        ini.Get(sec, "length_unit",
                ParametersLoadAccess::rgpot_options(params).length_unit);
    ParametersLoadAccess::rgpot_options(params).extensions_directory = ini.Get(
        sec, "extensions_directory",
        ParametersLoadAccess::rgpot_options(params).extensions_directory);
    ParametersLoadAccess::rgpot_options(params).check_consistency =
        ini.GetBoolean(
            sec, "check_consistency",
            ParametersLoadAccess::rgpot_options(params).check_consistency);
    ParametersLoadAccess::rgpot_options(params).uncertainty_threshold =
        ini.GetReal(
            sec, "uncertainty_threshold",
            ParametersLoadAccess::rgpot_options(params).uncertainty_threshold);
    ParametersLoadAccess::rgpot_options(params).torch_determinism_strict =
        ini.GetBoolean(sec, "torch_determinism_strict",
                       ParametersLoadAccess::rgpot_options(params)
                           .torch_determinism_strict);
    // XTB dlopen knobs (also accept [XTBPot] when backend=xtb)
    ParametersLoadAccess::rgpot_options(params).xtb_paramset = ini.Get(
        sec, "paramset",
        ini.Get(sec, "xtb_paramset",
                ParametersLoadAccess::rgpot_options(params).xtb_paramset));
    ParametersLoadAccess::rgpot_options(params).xtb_accuracy = ini.GetReal(
        sec, "accuracy",
        ini.GetReal(sec, "xtb_accuracy",
                    ParametersLoadAccess::rgpot_options(params).xtb_accuracy));
    ParametersLoadAccess::rgpot_options(params).xtb_electronic_temperature =
        ini.GetReal(sec, "electronic_temperature",
                    ini.GetReal(sec, "xtb_electronic_temperature",
                                ParametersLoadAccess::rgpot_options(params)
                                    .xtb_electronic_temperature));
    ParametersLoadAccess::rgpot_options(params).xtb_max_iterations =
        static_cast<int>(ini.GetInteger(
            sec, "max_iterations",
            ini.GetInteger(sec, "xtb_max_iterations",
                           ParametersLoadAccess::rgpot_options(params)
                               .xtb_max_iterations)));
    ParametersLoadAccess::rgpot_options(params).xtb_charge =
        ini.GetReal(sec, "xtb_charge",
                    static_cast<double>(
                        ParametersLoadAccess::rgpot_options(params).charge));
    ParametersLoadAccess::rgpot_options(params).xtb_uhf =
        static_cast<int>(ini.GetInteger(
            sec, "uhf",
            ini.GetInteger(
                sec, "xtb_uhf",
                ParametersLoadAccess::rgpot_options(params).xtb_uhf)));
    const std::string be =
        toLowerCase(ParametersLoadAccess::rgpot_options(params).backend);
    if (be == "xtb" || be == "xtbpot" || be == "gfn" || be == "gfnxtb") {
      ParametersLoadAccess::rgpot_options(params).xtb_paramset =
          ini.Get("XTBPot", "paramset",
                  ParametersLoadAccess::rgpot_options(params).xtb_paramset);
      ParametersLoadAccess::rgpot_options(params).xtb_accuracy =
          ini.GetReal("XTBPot", "accuracy",
                      ParametersLoadAccess::rgpot_options(params).xtb_accuracy);
      ParametersLoadAccess::rgpot_options(params).xtb_electronic_temperature =
          ini.GetReal("XTBPot", "electronic_temperature",
                      ParametersLoadAccess::rgpot_options(params)
                          .xtb_electronic_temperature);
      ParametersLoadAccess::rgpot_options(params).xtb_max_iterations =
          static_cast<int>(ini.GetInteger(
              "XTBPot", "max_iterations",
              ParametersLoadAccess::rgpot_options(params).xtb_max_iterations));
      ParametersLoadAccess::rgpot_options(params).xtb_uhf = static_cast<int>(
          ini.GetInteger("XTBPot", "uhf",
                         ParametersLoadAccess::rgpot_options(params).xtb_uhf));
      ParametersLoadAccess::rgpot_options(params).xtb_charge =
          ini.GetReal("XTBPot", "charge",
                      ParametersLoadAccess::rgpot_options(params).xtb_charge);
    }
  }

  // [Debug] //

  ParametersLoadAccess::debug_options(params).write_movies =
      ini.GetBoolean("Debug", "write_movies",
                     ParametersLoadAccess::debug_options(params).write_movies);
  ParametersLoadAccess::debug_options(params).write_movies_interval =
      ini.GetInteger(
          "Debug", "write_movies_interval",
          ParametersLoadAccess::debug_options(params).write_movies_interval);
  ParametersLoadAccess::debug_options(params).write_deprecated_outs =
      ini.GetBoolean(
          "Debug", "write_deprecated_outs",
          ParametersLoadAccess::debug_options(params).write_deprecated_outs);
  ParametersLoadAccess::debug_options(params).estimate_neb_eigenvalues =
      ini.GetBoolean(
          "Debug", "estimate_neb_eigenvalues",
          ParametersLoadAccess::debug_options(params).estimate_neb_eigenvalues);
  ParametersLoadAccess::debug_options(params).neb_mmf =
      toLowerCase(ini.Get("Debug", "neb_mmf_estimator",
                          ParametersLoadAccess::debug_options(params).neb_mmf));

  // [Structure Comparison] //

  ParametersLoadAccess::structure_comparison_options(params)
      .distance_difference =
      ini.GetReal("Structure Comparison", "distance_difference",
                  ParametersLoadAccess::structure_comparison_options(params)
                      .distance_difference);
  ParametersLoadAccess::structure_comparison_options(params).neighbor_cutoff =
      ini.GetReal("Structure Comparison", "neighbor_cutoff",
                  ParametersLoadAccess::structure_comparison_options(params)
                      .neighbor_cutoff);
  ParametersLoadAccess::structure_comparison_options(params).check_rotation =
      ini.GetBoolean("Structure Comparison", "check_rotation",
                     ParametersLoadAccess::structure_comparison_options(params)
                         .check_rotation);
  ParametersLoadAccess::structure_comparison_options(params).energy_difference =
      ini.GetReal("Structure Comparison", "energy_difference",
                  ParametersLoadAccess::structure_comparison_options(params)
                      .energy_difference);
  ParametersLoadAccess::structure_comparison_options(params)
      .indistinguishable_atoms =
      ini.GetBoolean("Structure Comparison", "indistinguishable_atoms",
                     ParametersLoadAccess::structure_comparison_options(params)
                         .indistinguishable_atoms);
  ParametersLoadAccess::structure_comparison_options(params)
      .remove_translation =
      ini.GetBoolean("Structure Comparison", "remove_translation",
                     ParametersLoadAccess::structure_comparison_options(params)
                         .remove_translation);

  // [Process Search] //

  ParametersLoadAccess::process_search_options(params).minimize_first =
      ini.GetBoolean(
          "Process Search", "minimize_first",
          ParametersLoadAccess::process_search_options(params).minimize_first);
  ParametersLoadAccess::process_search_options(params).minimization_offset =
      ini.GetReal("Process Search", "minimization_offset",
                  ParametersLoadAccess::process_search_options(params)
                      .minimization_offset);

  // [Optimizers] //
  auto inp_optMethod =
      magic_enum::enum_cast<OptType>(ini.Get("Optimizer", "opt_method", "none"),
                                     magic_enum::case_insensitive)
          .value_or(OptType::Unknown);
  if (inp_optMethod != OptType::None) {
    ParametersLoadAccess::optimizer_options(params).method = inp_optMethod;
  }

  ParametersLoadAccess::optimizer_options(params).convergence_metric =
      toLowerCase(ini.Get(
          "Optimizer", "convergence_metric",
          ParametersLoadAccess::optimizer_options(params).convergence_metric));
  if (auto label = eonc::helpers::convergenceMetricLabel(
          ParametersLoadAccess::optimizer_options(params).convergence_metric)) {
    ParametersLoadAccess::optimizer_options(params).convergence_metric_label =
        *label;
  } else {
    EONC_LOG_ERROR(
        "unknown convergence_metric {}",
        ParametersLoadAccess::optimizer_options(params).convergence_metric);
    error = 1;
  }

  if (ini.HasSection("Refine")) {
    ParametersLoadAccess::optimizer_options(params).refine.method =
        magic_enum::enum_cast<OptType>(ini.Get("Refine", "opt_method", ""),
                                       magic_enum::case_insensitive)
            .value_or(OptType::None);
    ParametersLoadAccess::optimizer_options(params).refine.threshold =
        ini.GetReal(
            "Refine", "threshold",
            ParametersLoadAccess::optimizer_options(params).refine.threshold);
  }

  ParametersLoadAccess::optimizer_options(params).converged_force = ini.GetReal(
      "Optimizer", "converged_force",
      ParametersLoadAccess::optimizer_options(params).converged_force);
  ParametersLoadAccess::optimizer_options(params).max_iterations =
      static_cast<size_t>(ini.GetInteger(
          "Optimizer", "max_iterations",
          ParametersLoadAccess::optimizer_options(params).max_iterations));
  ParametersLoadAccess::optimizer_options(params).max_move =
      ini.GetReal("Optimizer", "max_move",
                  ParametersLoadAccess::optimizer_options(params).max_move);
  // Handle each optimizer separately
  if (ini.HasSection("QuickMin")) {
    ParametersLoadAccess::optimizer_options(params).time_step_input =
        ini.GetReal(
            "QuickMin", "time_step",
            ParametersLoadAccess::optimizer_options(params).time_step_input);
    ParametersLoadAccess::optimizer_options(params).time_step =
        ParametersLoadAccess::optimizer_options(params).time_step_input /
        ParametersLoadAccess::constants(params).timeUnit;
    ParametersLoadAccess::optimizer_options(params).quickmin.steepest_descent =
        ini.GetBoolean("Optimizer", "qm_steepest_descent",
                       ParametersLoadAccess::optimizer_options(params)
                           .quickmin.steepest_descent);
  }
  if (ini.HasSection("FIRE")) {
    EONC_LOG_WARNING("Overwriting QuickMin timestep with Fire timestep!!");
    ParametersLoadAccess::optimizer_options(params).time_step_input =
        ini.GetReal(
            "FIRE", "time_step",
            ParametersLoadAccess::optimizer_options(params).time_step_input);
    ParametersLoadAccess::optimizer_options(params).time_step =
        ParametersLoadAccess::optimizer_options(params).time_step_input /
        ParametersLoadAccess::constants(params).timeUnit;
    ParametersLoadAccess::optimizer_options(params).max_time_step_input =
        ini.GetReal("FIRE", "time_step_max",
                    ParametersLoadAccess::optimizer_options(params)
                        .max_time_step_input);
    ParametersLoadAccess::optimizer_options(params).max_time_step =
        ParametersLoadAccess::optimizer_options(params).max_time_step_input /
        ParametersLoadAccess::constants(params).timeUnit;
  }
  // 2014 optbench INI puts lbfgs_* on [Optimizer]. Prefer [LBFGS] when present.
  {
    const char *lbfgs_sec = ini.HasSection("LBFGS") ? "LBFGS" : "Optimizer";
    ParametersLoadAccess::optimizer_options(params).lbfgs.memory =
        ini.GetInteger(
            lbfgs_sec, "lbfgs_memory",
            ParametersLoadAccess::optimizer_options(params).lbfgs.memory);
    ParametersLoadAccess::optimizer_options(params).lbfgs.inverse_curvature =
        ini.GetReal(lbfgs_sec, "lbfgs_inverse_curvature",
                    ParametersLoadAccess::optimizer_options(params)
                        .lbfgs.inverse_curvature);
    ParametersLoadAccess::optimizer_options(params)
        .lbfgs.max_inverse_curvature =
        ini.GetReal(lbfgs_sec, "lbfgs_max_inverse_curvature",
                    ParametersLoadAccess::optimizer_options(params)
                        .lbfgs.max_inverse_curvature);
    auto &lbfgs = ParametersLoadAccess::optimizer_options(params).lbfgs;
    lbfgs.auto_scale = ini.GetBoolean(lbfgs_sec, "lbfgs_auto_scale",
                                      lbfgs.auto_scale);
    lbfgs.angle_reset = ini.GetBoolean(lbfgs_sec, "lbfgs_angle_reset",
                                       lbfgs.angle_reset);
    lbfgs.distance_reset = ini.GetBoolean(lbfgs_sec, "lbfgs_distance_reset",
                                          lbfgs.distance_reset);
    lbfgs.curvature = toLowerCase(
        ini.Get(lbfgs_sec, "lbfgs_curvature", lbfgs.curvature));
    lbfgs.project_rigid = ini.GetBoolean(lbfgs_sec, "lbfgs_project_rigid",
                                         lbfgs.project_rigid);
    lbfgs.secant =
        toLowerCase(ini.Get(lbfgs_sec, "lbfgs_secant", lbfgs.secant));
    lbfgs.precon =
        toLowerCase(ini.Get(lbfgs_sec, "lbfgs_precon", lbfgs.precon));
    lbfgs.step = toLowerCase(ini.Get(lbfgs_sec, "lbfgs_step", lbfgs.step));
    lbfgs.h0 = toLowerCase(ini.Get(lbfgs_sec, "lbfgs_h0", lbfgs.h0));
    lbfgs.accept =
        toLowerCase(ini.Get(lbfgs_sec, "lbfgs_accept", lbfgs.accept));
    lbfgs.extra_updates = ini.GetInteger(lbfgs_sec, "lbfgs_extra_updates",
                                         lbfgs.extra_updates);
    lbfgs.cautious_eps =
        ini.GetReal(lbfgs_sec, "lbfgs_cautious_eps", lbfgs.cautious_eps);
    lbfgs.cautious_alpha =
        ini.GetReal(lbfgs_sec, "lbfgs_cautious_alpha", lbfgs.cautious_alpha);
    lbfgs.precon_A = ini.GetReal(lbfgs_sec, "lbfgs_precon_A", lbfgs.precon_A);
    lbfgs.precon_mu =
        ini.GetReal(lbfgs_sec, "lbfgs_precon_mu", lbfgs.precon_mu);
    lbfgs.precon_rcut =
        ini.GetReal(lbfgs_sec, "lbfgs_precon_rcut", lbfgs.precon_rcut);
  }
  {
    const char *xs = ini.HasSection("Xtsci") ? "Xtsci" : "Optimizer";
    params.optimizer_options.xtsci_method = toLowerCase(ini.Get(
        xs, "xtsci_method", params.optimizer_options.xtsci_method));
  }
  if (ini.HasSection("CG")) {
    ParametersLoadAccess::optimizer_options(params).cg.no_overshooting =
        ini.GetBoolean(
            "CG", "cg_no_overshooting",
            ParametersLoadAccess::optimizer_options(params).cg.no_overshooting);
    ParametersLoadAccess::optimizer_options(params).cg.knock_out_max_move =
        ini.GetBoolean("CG", "cg_knock_out_max_move",
                       ParametersLoadAccess::optimizer_options(params)
                           .cg.knock_out_max_move);
    ParametersLoadAccess::optimizer_options(params).cg.line_search =
        ini.GetBoolean(
            "CG", "cg_line_search",
            ParametersLoadAccess::optimizer_options(params).cg.line_search);
    ParametersLoadAccess::optimizer_options(params).cg.line_converged =
        ini.GetReal(
            "CG", "cg_line_converged",
            ParametersLoadAccess::optimizer_options(params).cg.line_converged);
    ParametersLoadAccess::optimizer_options(params).cg.max_iter_before_reset =
        ini.GetInteger("CG", "cg_max_iter_before_reset",
                       ParametersLoadAccess::optimizer_options(params)
                           .cg.max_iter_before_reset);
    ParametersLoadAccess::optimizer_options(params).cg.line_search_max_iter =
        ini.GetInteger("CG", "cg_max_iter_line_search",
                       ParametersLoadAccess::optimizer_options(params)
                           .cg.line_search_max_iter);
  }
  if (ini.HasSection("SD")) {
    ParametersLoadAccess::optimizer_options(params).sd.alpha =
        ini.GetReal("SD", "sd_alpha",
                    ParametersLoadAccess::optimizer_options(params).sd.alpha);
    ParametersLoadAccess::optimizer_options(params).sd.two_point =
        ini.GetBoolean(
            "SD", "sd_twopoint",
            ParametersLoadAccess::optimizer_options(params).sd.two_point);
  }

  // [Dimer] //

  ParametersLoadAccess::dimer_options(params).rotation_angle =
      ini.GetReal("Dimer", "finite_angle",
                  ParametersLoadAccess::dimer_options(params).rotation_angle);
  ParametersLoadAccess::dimer_options(params).improved =
      ini.GetBoolean("Dimer", "improved",
                     ParametersLoadAccess::dimer_options(params).improved);
  ParametersLoadAccess::dimer_options(params).converged_angle =
      ini.GetReal("Dimer", "converged_angle",
                  ParametersLoadAccess::dimer_options(params).converged_angle);
  ParametersLoadAccess::dimer_options(params).max_iterations = ini.GetInteger(
      "Dimer", "max_iterations",
      ParametersLoadAccess::dimer_options(params).max_iterations);
  if (auto dimerOpt = magic_enum::enum_cast<OptType>(
          ini.Get("Dimer", "opt_method", "cg"), magic_enum::case_insensitive);
      dimerOpt && *dimerOpt != OptType::Unknown && *dimerOpt != OptType::None) {
    ParametersLoadAccess::dimer_options(params).opt_method = *dimerOpt;
  }
  ParametersLoadAccess::dimer_options(params).rotations_min =
      ini.GetInteger("Dimer", "rotations_min",
                     ParametersLoadAccess::dimer_options(params).rotations_min);
  ParametersLoadAccess::dimer_options(params).rotations_max =
      ini.GetInteger("Dimer", "rotations_max",
                     ParametersLoadAccess::dimer_options(params).rotations_max);
  ParametersLoadAccess::dimer_options(params).torque_min =
      ini.GetReal("Dimer", "torque_min",
                  ParametersLoadAccess::dimer_options(params).torque_min);
  ParametersLoadAccess::dimer_options(params).torque_max =
      ini.GetReal("Dimer", "torque_max",
                  ParametersLoadAccess::dimer_options(params).torque_max);
  ParametersLoadAccess::dimer_options(params).remove_rotation = ini.GetBoolean(
      "Dimer", "remove_rotation",
      ParametersLoadAccess::dimer_options(params).remove_rotation);
  ParametersLoadAccess::dimer_options(params).lor_residual_tol =
      ini.GetReal("Dimer", "lor_residual_tol",
                  ParametersLoadAccess::dimer_options(params).lor_residual_tol);
  {
    const auto rotTok =
        toLowerCase(ini.Get("Dimer", "rotation_backend", "classical"));
    ParametersLoadAccess::dimer_options(params).rotation_backend =
        magic_enum::enum_cast<DimerRotationBackend>(
            rotTok, magic_enum::case_insensitive)
            .value_or(DimerRotationBackend::Classical);
  }

  // GP Surrogate Parameters
  ParametersLoadAccess::gp_surrogate_options(params).enabled =
      ini.GetBoolean("Surrogate", "use_surrogate", false);
  if (ParametersLoadAccess::gp_surrogate_options(params).enabled) {
    ParametersLoadAccess::gp_surrogate_options(params).sub_job =
        ParametersLoadAccess::main_options(params).job;
    ParametersLoadAccess::main_options(params).job = JobType::GP_Surrogate;
  }
  ParametersLoadAccess::gp_surrogate_options(params).uncertainty = ini.GetReal(
      "Surrogate", "gp_uncertainty",
      ParametersLoadAccess::gp_surrogate_options(params).uncertainty);
  if (ini.HasSection("Surrogate")) {
    ParametersLoadAccess::gp_surrogate_options(params).potential =
        magic_enum::enum_cast<PotType>(ini.Get("Surrogate", "potential", ""),
                                       magic_enum::case_insensitive)
            .value_or(PotType::UNKNOWN);
    if (params.gp_surrogate_options().potential != PotType::CatLearn) {
      throw std::runtime_error("We only support catlearn for GP right now");
    }
  }
  // [CatLearn]
  if (ini.HasSection("CatLearn")) {
    ParametersLoadAccess::catlearn_options(params).path =
        ini.Get("CatLearn", "catl_path", "");
    ParametersLoadAccess::catlearn_options(params).model =
        ini.Get("CatLearn", "model", "catl_model");
    ParametersLoadAccess::catlearn_options(params).prior =
        ini.Get("CatLearn", "prior", "catl_prior");
    ParametersLoadAccess::catlearn_options(params).use_deriv =
        ini.GetBoolean("CatLearn", "use_derivatives", "catl_deriv");
    ParametersLoadAccess::catlearn_options(params).use_fingerprint =
        ini.GetBoolean("CatLearn", "use_fingerprint", "catl_fingerprint");
    ParametersLoadAccess::catlearn_options(params).parallel = ini.GetBoolean(
        "CatLearn", "parallel_hyperparameter_opt", "catl_parallel");
  }
  // [ASE_ORCA]
  if (ini.HasSection("ASE_ORCA")) {
    ParametersLoadAccess::ase_orca_options(params).path =
        ini.Get("ASE_ORCA", "orca_path", "");
    ParametersLoadAccess::ase_orca_options(params).nproc =
        ini.Get("ASE_ORCA", "nproc", "1");
    ParametersLoadAccess::ase_orca_options(params).simpleinput =
        ini.Get("ASE_ORCA", "simpleinput", "");
    ParametersLoadAccess::ase_orca_options(params).charge =
        static_cast<int>(ini.GetInteger("ASE_ORCA", "charge", 0));
    ParametersLoadAccess::ase_orca_options(params).multiplicity =
        static_cast<int>(ini.GetInteger("ASE_ORCA", "multiplicity", 1));
  }
  // [ASE_NWCHEM]
  if (ini.HasSection("ASE_NWCHEM")) {
    ParametersLoadAccess::ase_nwchem_options(params).path =
        ini.Get("ASE_NWCHEM", "nwchem_path", "");
    ParametersLoadAccess::ase_nwchem_options(params).nproc =
        ini.Get("ASE_NWCHEM", "nproc", "1");
    ParametersLoadAccess::ase_nwchem_options(params).mpi_launcher =
        ini.Get("ASE_NWCHEM", "mpi_launcher",
                ParametersLoadAccess::ase_nwchem_options(params).mpi_launcher);
    ParametersLoadAccess::ase_nwchem_options(params).multiplicity =
        ini.Get("ASE_NWCHEM", "multiplicity", "");
    ParametersLoadAccess::ase_nwchem_options(params).scf_thresh =
        ini.GetReal("ASE_NWCHEM", "scf_thresh", 1e-5);
    ParametersLoadAccess::ase_nwchem_options(params).scf_maxiter =
        ini.GetInteger("ASE_NWCHEM", "scf_maxiter", 200);
    ParametersLoadAccess::ase_nwchem_options(params).basis =
        ini.Get("ASE_NWCHEM", "basis",
                ParametersLoadAccess::ase_nwchem_options(params).basis);
    ParametersLoadAccess::ase_nwchem_options(params).memory =
        ini.Get("ASE_NWCHEM", "memory",
                ParametersLoadAccess::ase_nwchem_options(params).memory);
  }
  // [Metatomic]
  if (ini.HasSection("Metatomic")) {
    ParametersLoadAccess::metatomic_options(params).model_path =
        ini.Get("Metatomic", "model_path", "");
    ParametersLoadAccess::metatomic_options(params).device =
        ini.Get("Metatomic", "device", "cpu");
    ParametersLoadAccess::metatomic_options(params).length_unit =
        ini.Get("Metatomic", "length_unit", "angstrom");
    ParametersLoadAccess::metatomic_options(params).extensions_directory =
        ini.Get("Metatomic", "extensions_directory", "");
    ParametersLoadAccess::metatomic_options(params).check_consistency =
        ini.GetBoolean("Metatomic", "check_consistency", false);
    ParametersLoadAccess::metatomic_options(params).uncertainty_threshold =
        ini.GetReal("Metatomic", "uncertainty_threshold", -1.0);
    ParametersLoadAccess::metatomic_options(params).energy_output =
        ini.Get("Metatomic", "energy_output",
                ParametersLoadAccess::metatomic_options(params).energy_output);
    ParametersLoadAccess::metatomic_options(params).energy_uncertainty_output =
        ini.Get("Metatomic", "energy_uncertainty_output",
                ParametersLoadAccess::metatomic_options(params)
                    .energy_uncertainty_output);
    ParametersLoadAccess::metatomic_options(params).force_output =
        ini.Get("Metatomic", "force_output",
                ParametersLoadAccess::metatomic_options(params).force_output);
    ParametersLoadAccess::metatomic_options(params).non_conservative =
        ini.GetBoolean("Metatomic", "non_conservative", false);
    ParametersLoadAccess::metatomic_options(params).random_rotation =
        ini.GetBoolean("Metatomic", "random_rotation", false);
    ParametersLoadAccess::metatomic_options(params).n_symmetry_rotations =
        static_cast<long>(
            ini.GetInteger("Metatomic", "n_symmetry_rotations", 0));
    ParametersLoadAccess::metatomic_options(params).deterministic =
        ini.GetBoolean("Metatomic", "deterministic", true);
    ParametersLoadAccess::metatomic_options(params).deterministic_strict =
        ini.GetBoolean("Metatomic", "deterministic_strict", false);
    auto &_variant = ParametersLoadAccess::metatomic_options(params).variant;
    _variant.base = ini.Get("Metatomic", "variant_base", "");
    _variant.energy = ini.Get("Metatomic", "variant_energy", "");
    _variant.energy_uncertainty =
        ini.Get("Metatomic", "variant_energy_uncertainty", "");
    _variant.force = ini.Get("Metatomic", "variant_force", "");
  }
  // [Serve]
  if (ini.HasSection("Serve")) {
    ParametersLoadAccess::serve_options(params).host = ini.Get(
        "Serve", "host", ParametersLoadAccess::serve_options(params).host);
    ParametersLoadAccess::serve_options(params).port =
        static_cast<uint16_t>(ini.GetInteger(
            "Serve", "port", ParametersLoadAccess::serve_options(params).port));
    ParametersLoadAccess::serve_options(params).replicas = static_cast<size_t>(
        ini.GetInteger("Serve", "replicas",
                       ParametersLoadAccess::serve_options(params).replicas));
    ParametersLoadAccess::serve_options(params).gateway_port =
        static_cast<uint16_t>(ini.GetInteger(
            "Serve", "gateway_port",
            ParametersLoadAccess::serve_options(params).gateway_port));
    ParametersLoadAccess::serve_options(params).endpoints =
        ini.Get("Serve", "endpoints",
                ParametersLoadAccess::serve_options(params).endpoints);
  }

  // GP_NEB only
  ParametersLoadAccess::gp_surrogate_options(params).linear_path_always =
      ini.GetBoolean("Surrogate", "gp_linear_path_always",
                     ParametersLoadAccess::gp_surrogate_options(params)
                         .linear_path_always);
  // [Lanczos] //

  ParametersLoadAccess::lanczos_options(params).tolerance =
      ini.GetReal("Lanczos", "tolerance",
                  ParametersLoadAccess::lanczos_options(params).tolerance);
  ParametersLoadAccess::lanczos_options(params).max_iterations = ini.GetInteger(
      "Lanczos", "max_iterations",
      ParametersLoadAccess::lanczos_options(params).max_iterations);
  ParametersLoadAccess::lanczos_options(params).quit_early =
      ini.GetBoolean("Lanczos", "quit_early",
                     ParametersLoadAccess::lanczos_options(params).quit_early);
  if (ini.HasValue("Lanczos", "phva_atoms")) {
    ParametersLoadAccess::lanczos_options(params).phva_atoms =
        toLowerCase(ini.Get("Lanczos", "phva_atoms", "All"));
  }

  // [Davidson] //
  ParametersLoadAccess::davidson_options(params).tolerance =
      ini.GetReal("Davidson", "tolerance",
                  ParametersLoadAccess::davidson_options(params).tolerance);
  ParametersLoadAccess::davidson_options(params).max_iterations =
      ini.GetInteger(
          "Davidson", "max_iterations",
          ParametersLoadAccess::davidson_options(params).max_iterations);
  ParametersLoadAccess::davidson_options(params).diagonal_preconditioner =
      ini.GetBoolean("Davidson", "diagonal_preconditioner",
                     ParametersLoadAccess::davidson_options(params)
                         .diagonal_preconditioner);
  if (ini.HasValue("Davidson", "phva_atoms")) {
    ParametersLoadAccess::davidson_options(params).phva_atoms =
        toLowerCase(ini.Get("Davidson", "phva_atoms", "All"));
  }

  // [ARTn] //
  ParametersLoadAccess::artn_options(params).push_step_size =
      ini.GetReal("ARTn", "push_step_size",
                  ParametersLoadAccess::artn_options(params).push_step_size);
  ParametersLoadAccess::artn_options(params).force_threshold =
      ini.GetReal("ARTn", "force_threshold",
                  ParametersLoadAccess::artn_options(params).force_threshold);
  ParametersLoadAccess::artn_options(params).max_iterations =
      ini.GetInteger("ARTn", "max_iterations",
                     ParametersLoadAccess::artn_options(params).max_iterations);
  ParametersLoadAccess::artn_options(params).ninit = ini.GetInteger(
      "ARTn", "ninit", ParametersLoadAccess::artn_options(params).ninit);
  ParametersLoadAccess::artn_options(params).nperp_limitation =
      ini.Get("ARTn", "nperp_limitation",
              ParametersLoadAccess::artn_options(params).nperp_limitation);
  ParametersLoadAccess::artn_options(params).lanczos_min_size = ini.GetInteger(
      "ARTn", "lanczos_min_size",
      ParametersLoadAccess::artn_options(params).lanczos_min_size);
  ParametersLoadAccess::artn_options(params).nsmooth = ini.GetInteger(
      "ARTn", "nsmooth", ParametersLoadAccess::artn_options(params).nsmooth);
  ParametersLoadAccess::artn_options(params).filin = ini.Get(
      "ARTn", "filin", ParametersLoadAccess::artn_options(params).filin);

  // [IRA] //
  ParametersLoadAccess::ira_options(params).distance_threshold =
      ini.GetReal("IRA", "distance_threshold",
                  ParametersLoadAccess::ira_options(params).distance_threshold);
  ParametersLoadAccess::ira_options(params).symmetry_threshold =
      ini.GetReal("IRA", "symmetry_threshold",
                  ParametersLoadAccess::ira_options(params).symmetry_threshold);
  ParametersLoadAccess::ira_options(params).use_pbc = ini.GetBoolean(
      "IRA", "use_pbc", ParametersLoadAccess::ira_options(params).use_pbc);

  // [GPR Dimer] //
  ParametersLoadAccess::gpr_dimer_options(params).rotation_angle = ini.GetReal(
      "GPR Dimer", "finite_angle",
      ParametersLoadAccess::gpr_dimer_options(params).rotation_angle);
  ParametersLoadAccess::gpr_dimer_options(params).converged_angle = ini.GetReal(
      "GPR Dimer", "converged_angle",
      ParametersLoadAccess::gpr_dimer_options(params).converged_angle);
  ParametersLoadAccess::gpr_dimer_options(params).relax_conv_angle =
      ini.GetReal(
          "GPR Dimer", "relaxation_converged_angle",
          ParametersLoadAccess::gpr_dimer_options(params).relax_conv_angle);
  ParametersLoadAccess::gpr_dimer_options(params).init_rotations_max =
      static_cast<long>(ini.GetInteger(
          "GPR Dimer", "max_initial_rotation_iterations",
          ParametersLoadAccess::gpr_dimer_options(params).init_rotations_max));
  ParametersLoadAccess::gpr_dimer_options(params).relax_rotations_max =
      static_cast<long>(ini.GetInteger(
          "GPR Dimer", "max_relaxation_rotation_iterations",
          ParametersLoadAccess::gpr_dimer_options(params).relax_rotations_max));
  ParametersLoadAccess::gpr_dimer_options(params).divisor_t_dimer_gp =
      static_cast<long>(ini.GetInteger(
          "GPR Dimer", "divisor_t_dimer",
          ParametersLoadAccess::gpr_dimer_options(params).divisor_t_dimer_gp));
  ParametersLoadAccess::gpr_dimer_options(params)
      .max_outer_iterations = static_cast<long>(ini.GetInteger(
      "GPR Dimer", "max_outer_iterations",
      ParametersLoadAccess::gpr_dimer_options(params).max_outer_iterations));
  ParametersLoadAccess::gpr_dimer_options(params)
      .max_inner_iterations = static_cast<long>(ini.GetInteger(
      "GPR Dimer", "max_inner_iterations",
      ParametersLoadAccess::gpr_dimer_options(params).max_inner_iterations));
  ParametersLoadAccess::gpr_dimer_options(params).midpoint_max_disp =
      ini.GetReal(
          "GPR Dimer", "max_midpoint_displacement",
          ParametersLoadAccess::gpr_dimer_options(params).midpoint_max_disp);
  ParametersLoadAccess::gpr_dimer_options(params).rot_opt_method =
      ini.Get("GPR Dimer", "rotation_opt_method",
              ParametersLoadAccess::gpr_dimer_options(params).rot_opt_method);
  ParametersLoadAccess::gpr_dimer_options(params).trans_opt_method =
      ini.Get("GPR Dimer", "translation_opt_method",
              ParametersLoadAccess::gpr_dimer_options(params).trans_opt_method);
  ParametersLoadAccess::gpr_dimer_options(params).active_radius = ini.GetReal(
      "GPR Dimer", "active_radius",
      ParametersLoadAccess::gpr_dimer_options(params).active_radius);
  ParametersLoadAccess::gpr_dimer_options(params).dimer_sep =
      ini.GetReal("GPR Dimer", "dimer_separation",
                  ParametersLoadAccess::gpr_dimer_options(params).dimer_sep);
  ParametersLoadAccess::gpr_dimer_options(params).conv_step =
      ini.GetReal("GPR Dimer", "convex_region_step_size",
                  ParametersLoadAccess::gpr_dimer_options(params).conv_step);
  ParametersLoadAccess::gpr_dimer_options(params).max_step =
      ini.GetReal("GPR Dimer", "max_step_size",
                  ParametersLoadAccess::gpr_dimer_options(params).max_step);
  ParametersLoadAccess::gpr_dimer_options(params).ratio_at_limit = ini.GetReal(
      "GPR Dimer", "ratio_at_limit",
      ParametersLoadAccess::gpr_dimer_options(params).ratio_at_limit);
  ParametersLoadAccess::gpr_dimer_options(params).init_rot_gp = ini.GetBoolean(
      "GPR Dimer", "nogp_initial_rotations",
      ParametersLoadAccess::gpr_dimer_options(params).init_rot_gp);
  ParametersLoadAccess::gpr_dimer_options(params).init_trans_gp =
      ini.GetBoolean(
          "GPR Dimer", "nogp_init_translations",
          ParametersLoadAccess::gpr_dimer_options(params).init_trans_gp);
  ParametersLoadAccess::gpr_dimer_options(params).many_iterations =
      ini.GetBoolean(
          "GPR Dimer", "has_many_iterations",
          ParametersLoadAccess::gpr_dimer_options(params).many_iterations);
  // GPR Params
  ParametersLoadAccess::gpr_dimer_options(params).gpr_params.hyper_opt_method =
      ini.Get("GPR Dimer", "hyperparameter_opt_method",
              ParametersLoadAccess::gpr_dimer_options(params)
                  .gpr_params.hyper_opt_method);
  ParametersLoadAccess::gpr_dimer_options(params).gpr_params.sigma2 =
      ini.GetReal(
          "GPR Dimer", "gpr_variance",
          ParametersLoadAccess::gpr_dimer_options(params).gpr_params.sigma2);
  ParametersLoadAccess::gpr_dimer_options(params).gpr_params.jitter_sigma2 =
      ini.GetReal("GPR Dimer", "gpr_jitter_variance",
                  ParametersLoadAccess::gpr_dimer_options(params)
                      .gpr_params.jitter_sigma2);
  ParametersLoadAccess::gpr_dimer_options(params).gpr_params.noise_sigma2 =
      ini.GetReal("GPR Dimer", "gpr_noise_variance",
                  ParametersLoadAccess::gpr_dimer_options(params)
                      .gpr_params.noise_sigma2);
  ParametersLoadAccess::gpr_dimer_options(params).gpr_params.prior_mu =
      ini.GetReal(
          "GPR Dimer", "prior_mean",
          ParametersLoadAccess::gpr_dimer_options(params).gpr_params.prior_mu);
  ParametersLoadAccess::gpr_dimer_options(params).gpr_params.prior_sigma2 =
      ini.GetReal("GPR Dimer", "prior_variance",
                  ParametersLoadAccess::gpr_dimer_options(params)
                      .gpr_params.prior_sigma2);
  ParametersLoadAccess::gpr_dimer_options(params).gpr_params.prior_nu =
      ini.GetReal(
          "GPR Dimer", "prior_degrees_of_freedom",
          ParametersLoadAccess::gpr_dimer_options(params).gpr_params.prior_nu);
  // GPR Optimization Parameters
  ParametersLoadAccess::gpr_dimer_options(params).opt_params.check_derivatives =
      ini.GetBoolean("GPR Dimer", "check_derivatives",
                     ParametersLoadAccess::gpr_dimer_options(params)
                         .opt_params.check_derivatives);
  ParametersLoadAccess::gpr_dimer_options(params).opt_params.max_iterations =
      static_cast<int>(
          ini.GetInteger("GPR Dimer", "opt_max_iterations",
                         ParametersLoadAccess::gpr_dimer_options(params)
                             .opt_params.max_iterations));
  ParametersLoadAccess::gpr_dimer_options(params).opt_params.tol_func =
      ini.GetReal(
          "GPR Dimer", "opt_tol_func",
          ParametersLoadAccess::gpr_dimer_options(params).opt_params.tol_func);
  ParametersLoadAccess::gpr_dimer_options(params).opt_params.tol_sol =
      ini.GetReal(
          "GPR Dimer", "opt_tol_sol",
          ParametersLoadAccess::gpr_dimer_options(params).opt_params.tol_sol);
  ParametersLoadAccess::gpr_dimer_options(params).opt_params.lambda_limit =
      ini.GetReal("GPR Dimer", "opt_lambda_limit",
                  ParametersLoadAccess::gpr_dimer_options(params)
                      .opt_params.lambda_limit);
  ParametersLoadAccess::gpr_dimer_options(params).opt_params.lambda_init =
      ini.GetReal("GPR Dimer", "opt_lambda_init",
                  ParametersLoadAccess::gpr_dimer_options(params)
                      .opt_params.lambda_init);
  // GPR Debugging Parameters
  ParametersLoadAccess::gpr_dimer_options(params).debug_params.report_level =
      static_cast<int>(
          ini.GetInteger("GPR Dimer", "report_level",
                         ParametersLoadAccess::gpr_dimer_options(params)
                             .debug_params.report_level));
  ParametersLoadAccess::gpr_dimer_options(params).debug_params.debug_level =
      static_cast<int>(
          ini.GetInteger("GPR Dimer", "debug_level",
                         ParametersLoadAccess::gpr_dimer_options(params)
                             .debug_params.debug_level));
  ParametersLoadAccess::gpr_dimer_options(params).debug_params.out_dir =
      ini.Get(
          "GPR Dimer", "debug_output_directory",
          ParametersLoadAccess::gpr_dimer_options(params).debug_params.out_dir);
  ParametersLoadAccess::gpr_dimer_options(params).debug_params.pos_file =
      ini.Get("GPR Dimer", "debug_position_basename",
              ParametersLoadAccess::gpr_dimer_options(params)
                  .debug_params.pos_file);
  ParametersLoadAccess::gpr_dimer_options(params).debug_params.energy_file =
      ini.Get("GPR Dimer", "debug_energy_basename",
              ParametersLoadAccess::gpr_dimer_options(params)
                  .debug_params.energy_file);
  ParametersLoadAccess::gpr_dimer_options(params).debug_params.grad_file =
      ini.Get("GPR Dimer", "debug_gradient_basename",
              ParametersLoadAccess::gpr_dimer_options(params)
                  .debug_params.grad_file);
  ParametersLoadAccess::gpr_dimer_options(params)
      .debug_params.offset_mid_point =
      ini.GetReal("GPR Dimer", "debug_midpoint_offset",
                  ParametersLoadAccess::gpr_dimer_options(params)
                      .debug_params.offset_mid_point);
  ParametersLoadAccess::gpr_dimer_options(params).debug_params.dy = ini.GetReal(
      "GPR Dimer", "debug_y_step",
      ParametersLoadAccess::gpr_dimer_options(params).debug_params.dy);
  ParametersLoadAccess::gpr_dimer_options(params).debug_params.dz = ini.GetReal(
      "GPR Dimer", "debug_z_step",
      ParametersLoadAccess::gpr_dimer_options(params).debug_params.dz);
  // GPR Prune
  ParametersLoadAccess::gpr_dimer_options(params).prune_params.use_prune =
      ini.GetBoolean("GPR Dimer", "use_prune",
                     ParametersLoadAccess::gpr_dimer_options(params)
                         .prune_params.use_prune);
  ParametersLoadAccess::gpr_dimer_options(params).prune_params.begin =
      static_cast<int>(ini.GetInteger(
          "GPR Dimer", "start_prune_at",
          ParametersLoadAccess::gpr_dimer_options(params).prune_params.begin));
  ParametersLoadAccess::gpr_dimer_options(params).prune_params.n_vals =
      static_cast<int>(ini.GetInteger(
          "GPR Dimer", "nprune_vals",
          ParametersLoadAccess::gpr_dimer_options(params).prune_params.n_vals));
  ParametersLoadAccess::gpr_dimer_options(params).prune_params.threshold =
      ini.GetReal("GPR Dimer", "prune_threshold",
                  ParametersLoadAccess::gpr_dimer_options(params)
                      .prune_params.threshold);

  // [Prefactor] //

  ParametersLoadAccess::prefactor_options(params).default_value = ini.GetReal(
      "Prefactor", "default_value",
      ParametersLoadAccess::prefactor_options(params).default_value);
  ParametersLoadAccess::prefactor_options(params).max_value =
      ini.GetReal("Prefactor", "max_value",
                  ParametersLoadAccess::prefactor_options(params).max_value);
  ParametersLoadAccess::prefactor_options(params).min_value =
      ini.GetReal("Prefactor", "min_value",
                  ParametersLoadAccess::prefactor_options(params).min_value);
  ParametersLoadAccess::prefactor_options(params).within_radius = ini.GetReal(
      "Prefactor", "within_radius",
      ParametersLoadAccess::prefactor_options(params).within_radius);
  ParametersLoadAccess::prefactor_options(params).min_displacement =
      ini.GetReal(
          "Prefactor", "min_displacement",
          ParametersLoadAccess::prefactor_options(params).min_displacement);
  ParametersLoadAccess::prefactor_options(params).rate = toLowerCase(
      ini.Get("Prefactor", "rate_estimation",
              ParametersLoadAccess::prefactor_options(params).rate));
  ParametersLoadAccess::prefactor_options(params).configuration = toLowerCase(
      ini.Get("Prefactor", "configuration",
              ParametersLoadAccess::prefactor_options(params).configuration));
  ParametersLoadAccess::prefactor_options(params).all_free_atoms =
      ini.GetBoolean(
          "Prefactor", "all_free_atoms",
          ParametersLoadAccess::prefactor_options(params).all_free_atoms);
  ParametersLoadAccess::prefactor_options(params).filter_scheme = toLowerCase(
      ini.Get("Prefactor", "filter_scheme",
              ParametersLoadAccess::prefactor_options(params).filter_scheme));
  ParametersLoadAccess::prefactor_options(params).filter_fraction = ini.GetReal(
      "Prefactor", "filter_fraction",
      ParametersLoadAccess::prefactor_options(params).filter_fraction);

  // [Hessian] //
  // Prefer phva_atoms; accept legacy atom_list when phva_atoms is absent.
  if (ini.HasValue("Hessian", "phva_atoms")) {
    ParametersLoadAccess::hessian_options(params).phva_atoms =
        toLowerCase(ini.Get("Hessian", "phva_atoms", "All"));
  } else if (ini.HasValue("Hessian", "atom_list")) {
    ParametersLoadAccess::hessian_options(params).phva_atoms =
        toLowerCase(ini.Get("Hessian", "atom_list", "All"));
  }
  ParametersLoadAccess::hessian_options(params).zero_freq_value = ini.GetReal(
      "Hessian", "zero_freq_value",
      ParametersLoadAccess::hessian_options(params).zero_freq_value);
  ParametersLoadAccess::hessian_options(params).fd_scheme = toLowerCase(
      ini.Get("Hessian", "fd_scheme",
              ParametersLoadAccess::hessian_options(params).fd_scheme));
  ParametersLoadAccess::hessian_options(params).resume =
      ini.GetBoolean("Hessian", "resume",
                     ParametersLoadAccess::hessian_options(params).resume);
  ParametersLoadAccess::hessian_options(params).checkpoint_path =
      ini.Get("Hessian", "checkpoint_path",
              ParametersLoadAccess::hessian_options(params).checkpoint_path);

  // [Nudged Elastic Band] //
  const std::string neb_section = "Nudged Elastic Band";

  ParametersLoadAccess::neb_options(params).image_count =
      ini.GetInteger(neb_section, "images",
                     ParametersLoadAccess::neb_options(params).image_count);
  ParametersLoadAccess::neb_options(params).max_iterations = ini.GetInteger(
      neb_section, "max_iterations",
      ParametersLoadAccess::optimizer_options(params).max_iterations);
  ParametersLoadAccess::neb_options(params).force_tolerance = ini.GetReal(
      neb_section, "converged_force",
      ParametersLoadAccess::optimizer_options(params).converged_force);
  auto neb_optMethod =
      magic_enum::enum_cast<OptType>(ini.Get(neb_section, "opt_method", "none"),
                                     magic_enum::case_insensitive)
          .value_or(OptType::Unknown);
  if (neb_optMethod != OptType::None) {
    ParametersLoadAccess::neb_options(params).opt_method = neb_optMethod;
  }
  ParametersLoadAccess::neb_options(params).mmf_peaks.enabled = ini.GetBoolean(
      neb_section, "setup_mmf_peaks",
      ParametersLoadAccess::neb_options(params).mmf_peaks.enabled);
  ParametersLoadAccess::neb_options(params).mmf_peaks.tolerance = ini.GetReal(
      neb_section, "mmf_peak_tolerance",
      ParametersLoadAccess::neb_options(params).mmf_peaks.tolerance);
  ParametersLoadAccess::neb_options(params).match_endpoints =
      ini.GetBoolean(neb_section, "match_endpoints",
                     ParametersLoadAccess::neb_options(params).match_endpoints);
  ParametersLoadAccess::neb_options(params).match_method = toLowerCase(
      ini.Get(neb_section, "match_method",
              ParametersLoadAccess::neb_options(params).match_method));

  ParametersLoadAccess::neb_options(params).spring.constant =
      ini.GetReal(neb_section, "spring",
                  ParametersLoadAccess::neb_options(params).spring.constant);
  ParametersLoadAccess::neb_options(params).spring.use_elastic_band =
      ini.GetBoolean(
          neb_section, "elastic_band",
          ParametersLoadAccess::neb_options(params).spring.use_elastic_band);
  ParametersLoadAccess::neb_options(params).spring.doubly_nudged =
      ini.GetBoolean(
          neb_section, "doubly_nudged",
          ParametersLoadAccess::neb_options(params).spring.doubly_nudged);
  ParametersLoadAccess::neb_options(params).spring.use_switching =
      ini.GetBoolean(
          neb_section, "doubly_nudged_switching",
          ParametersLoadAccess::neb_options(params).spring.use_switching);

  ParametersLoadAccess::neb_options(params).spring.weighting.enabled =
      ini.GetBoolean(
          neb_section, "energy_weighted",
          ParametersLoadAccess::neb_options(params).spring.weighting.enabled);
  ParametersLoadAccess::neb_options(params).spring.weighting.trigger =
      ini.GetReal(
          neb_section, "ew_trigger",
          ParametersLoadAccess::neb_options(params).spring.weighting.trigger);
  ParametersLoadAccess::neb_options(params).spring.weighting.k_min =
      ini.GetReal(
          neb_section, "ew_ksp_min",
          ParametersLoadAccess::neb_options(params).spring.weighting.k_min);
  ParametersLoadAccess::neb_options(params).spring.weighting.k_max =
      ini.GetReal(
          neb_section, "ew_ksp_max",
          ParametersLoadAccess::neb_options(params).spring.weighting.k_max);

  ParametersLoadAccess::neb_options(params).spring.om.enabled = ini.GetBoolean(
      neb_section, "onsager_machlup",
      ParametersLoadAccess::neb_options(params).spring.om.enabled);
  ParametersLoadAccess::neb_options(params).spring.om.optimize_k =
      ini.GetBoolean(
          neb_section, "om_optimize_k",
          ParametersLoadAccess::neb_options(params).spring.om.optimize_k);
  ParametersLoadAccess::neb_options(params).spring.om.k_scale =
      ini.GetReal(neb_section, "om_k_scale",
                  ParametersLoadAccess::neb_options(params).spring.om.k_scale);
  ParametersLoadAccess::neb_options(params).spring.om.k_min =
      ini.GetReal(neb_section, "om_k_min",
                  ParametersLoadAccess::neb_options(params).spring.om.k_min);
  ParametersLoadAccess::neb_options(params).spring.om.k_max =
      ini.GetReal(neb_section, "om_k_max",
                  ParametersLoadAccess::neb_options(params).spring.om.k_max);

  ParametersLoadAccess::neb_options(params).climbing_image.enabled =
      ini.GetBoolean(
          neb_section, "climbing_image_method",
          ParametersLoadAccess::neb_options(params).climbing_image.enabled);
  ParametersLoadAccess::neb_options(params).climbing_image.converged_only =
      ini.GetBoolean(neb_section, "climbing_image_converged_only",
                     ParametersLoadAccess::neb_options(params)
                         .climbing_image.converged_only);
  ParametersLoadAccess::neb_options(params).climbing_image.band_slack =
      ini.GetReal(
          neb_section, "climbing_image_band_slack",
          ParametersLoadAccess::neb_options(params).climbing_image.band_slack);
  ParametersLoadAccess::neb_options(params).climbing_image.use_old_tangent =
      ini.GetBoolean(neb_section, "old_tangent",
                     ParametersLoadAccess::neb_options(params)
                         .climbing_image.use_old_tangent);
  ParametersLoadAccess::neb_options(params).climbing_image.trigger_force =
      ini.GetReal(neb_section, "ci_after",
                  ParametersLoadAccess::neb_options(params)
                      .climbing_image.trigger_force);
  ParametersLoadAccess::neb_options(params).climbing_image.trigger_factor =
      ini.GetReal(neb_section, "ci_after_rel",
                  ParametersLoadAccess::neb_options(params)
                      .climbing_image.trigger_factor);

  auto &oci = ParametersLoadAccess::neb_options(params).climbing_image.ocineb;
  oci.use_mmf = ini.GetBoolean(neb_section, "ci_mmf", oci.use_mmf);
  oci.trigger_force =
      ini.GetReal(neb_section, "ci_mmf_after", oci.trigger_force);
  oci.trigger_factor =
      ini.GetReal(neb_section, "ci_mmf_after_rel", oci.trigger_factor);
  oci.max_steps = ini.GetInteger(neb_section, "ci_mmf_nsteps", oci.max_steps);
  oci.ci_stability_count = ini.GetInteger(
      neb_section, "ci_mmf_ci_stability_count", oci.ci_stability_count);
  oci.angle_tol = ini.GetReal(neb_section, "ci_mmf_angle", oci.angle_tol);

  auto &init = ParametersLoadAccess::neb_options(params).initialization;
  init.method =
      magic_enum::enum_cast<NEBInit>(ini.Get(neb_section, "initializer", ""),
                                     magic_enum::case_insensitive)
          .value_or(NEBInit::LINEAR);
  init.input_path = ini.Get(neb_section, "initial_path_in", init.input_path);
  init.max_iterations =
      ini.GetInteger(neb_section, "init_max_iterations", init.max_iterations);
  init.nsteps = ini.GetInteger(neb_section, "init_nsteps", init.nsteps);
  init.max_move = ini.GetReal(neb_section, "init_max_move", init.max_move);
  init.force_tolerance =
      ini.GetReal(neb_section, "init_force_threshold", init.force_tolerance);
  init.sidpp_alpha =
      ini.GetReal(neb_section, "sidpp_growth_alpha", init.sidpp_alpha);
  init.sidpp_frontier_tol =
      ini.GetReal(neb_section, "sidpp_frontier_tol", init.sidpp_frontier_tol);
  init.sidpp_reparam =
      ini.GetBoolean(neb_section, "sidpp_reparameterize", init.sidpp_reparam);
  init.sidpp_ideal_ksp =
      ini.GetBoolean(neb_section, "sidpp_ideal_ksp", init.sidpp_ideal_ksp);
  auto neb_ipath_optMethod =
      magic_enum::enum_cast<OptType>(
          ini.Get(neb_section, "ipath_opt_method", "none"),
          magic_enum::case_insensitive)
          .value_or(OptType::Unknown);
  if (neb_ipath_optMethod != OptType::None) {
    ParametersLoadAccess::neb_options(params).initialization.opt_method =
        neb_ipath_optMethod;
  }
  init.oversampling =
      ini.GetBoolean(neb_section, "oversampling", init.oversampling);
  init.oversampling_factor = ini.GetInteger(neb_section, "oversampling_factor",
                                            init.oversampling_factor);

  ParametersLoadAccess::neb_options(params).endpoints.minimize = ini.GetBoolean(
      neb_section, "minimize_endpoints",
      ParametersLoadAccess::neb_options(params).endpoints.minimize);
  ParametersLoadAccess::neb_options(params).endpoints.use_path_file =
      ini.GetBoolean(
          neb_section, "minimize_endpoints_for_ipath",
          ParametersLoadAccess::neb_options(params).endpoints.use_path_file);

  // [Dynamics] //

  ParametersLoadAccess::dynamics_options(params).time_step_input = ini.GetReal(
      "Dynamics", "time_step",
      ParametersLoadAccess::dynamics_options(params).time_step_input);
  ParametersLoadAccess::dynamics_options(params).time_step =
      ParametersLoadAccess::dynamics_options(params).time_step_input /
      ParametersLoadAccess::constants(params).timeUnit;
  ParametersLoadAccess::dynamics_options(params).time_input =
      ini.GetReal("Dynamics", "time",
                  ParametersLoadAccess::dynamics_options(params).time_input);
  ParametersLoadAccess::dynamics_options(params).time =
      ParametersLoadAccess::dynamics_options(params).time_input /
      ParametersLoadAccess::constants(params).timeUnit;
  if (ParametersLoadAccess::dynamics_options(params).time_step > 0.0) {
    ParametersLoadAccess::dynamics_options(params).steps =
        static_cast<long>(std::floor(
            ParametersLoadAccess::dynamics_options(params).time /
                ParametersLoadAccess::dynamics_options(params).time_step +
            0.5));
  } else {
    ParametersLoadAccess::dynamics_options(params).steps = 0;
  }
  ParametersLoadAccess::thermostat_options(params).kind =
      toLowerCase(ini.Get("Dynamics", "thermostat", "andersen"));
  ParametersLoadAccess::thermostat_options(params).andersen_alpha = ini.GetReal(
      "Dynamics", "andersen_alpha",
      ParametersLoadAccess::thermostat_options(params).andersen_alpha);
  ParametersLoadAccess::thermostat_options(params).andersen_tcol_input =
      ini.GetReal(
          "Dynamics", "andersen_collision_period",
          ParametersLoadAccess::thermostat_options(params).andersen_tcol_input);
  ParametersLoadAccess::thermostat_options(params).andersen_tcol =
      ParametersLoadAccess::thermostat_options(params).andersen_tcol_input /
      ParametersLoadAccess::constants(params).timeUnit;
  ParametersLoadAccess::thermostat_options(params).nose_mass =
      ini.GetReal("Dynamics", "nose_mass",
                  ParametersLoadAccess::thermostat_options(params).nose_mass);
  ParametersLoadAccess::thermostat_options(params).langevin_friction_input =
      ini.GetReal("Dynamics", "langevin_friction",
                  ParametersLoadAccess::thermostat_options(params)
                      .langevin_friction_input);
  ParametersLoadAccess::thermostat_options(params).langevin_friction =
      ParametersLoadAccess::thermostat_options(params).langevin_friction_input *
      ParametersLoadAccess::constants(params).timeUnit;

  // [Parallel Replica]

  ParametersLoadAccess::parallel_replica_options(params).auto_stop =
      ini.GetBoolean(
          "Parallel Replica", "stop_after_transition",
          ParametersLoadAccess::parallel_replica_options(params).auto_stop);
  ParametersLoadAccess::parallel_replica_options(params).refine_transition =
      ini.GetBoolean("Parallel Replica", "refine_transition",
                     ParametersLoadAccess::parallel_replica_options(params)
                         .refine_transition);
  ParametersLoadAccess::parallel_replica_options(params).dephase_loop_stop =
      ini.GetBoolean("Parallel Replica", "dephase_loop_stop",
                     ParametersLoadAccess::parallel_replica_options(params)
                         .dephase_loop_stop);
  ParametersLoadAccess::parallel_replica_options(params).dephase_time_input =
      ini.GetReal("Parallel Replica", "dephase_time",
                  ParametersLoadAccess::parallel_replica_options(params)
                      .dephase_time_input);
  ParametersLoadAccess::parallel_replica_options(params).dephase_time =
      ParametersLoadAccess::parallel_replica_options(params)
          .dephase_time_input /
      ParametersLoadAccess::constants(params).timeUnit;
  ParametersLoadAccess::parallel_replica_options(params).dephase_loop_max =
      ini.GetInteger("Parallel Replica", "dephase_loop_max",
                     ParametersLoadAccess::parallel_replica_options(params)
                         .dephase_loop_max);
  ParametersLoadAccess::parallel_replica_options(params)
      .state_check_interval_input =
      ini.GetReal("Parallel Replica", "state_check_interval",
                  ParametersLoadAccess::parallel_replica_options(params)
                      .state_check_interval_input);
  ParametersLoadAccess::parallel_replica_options(params).state_check_interval =
      ParametersLoadAccess::parallel_replica_options(params)
          .state_check_interval_input /
      ParametersLoadAccess::constants(params).timeUnit;
  ParametersLoadAccess::parallel_replica_options(params).record_interval_input =
      ini.GetReal("Parallel Replica", "state_save_interval",
                  0.1 * ParametersLoadAccess::parallel_replica_options(params)
                            .state_check_interval_input);
  ParametersLoadAccess::parallel_replica_options(params).record_interval =
      ParametersLoadAccess::parallel_replica_options(params)
          .record_interval_input /
      ParametersLoadAccess::constants(params).timeUnit;
  ParametersLoadAccess::parallel_replica_options(params).corr_time_input =
      ini.GetReal("Parallel Replica", "post_transition_time",
                  ParametersLoadAccess::parallel_replica_options(params)
                      .corr_time_input);
  ParametersLoadAccess::parallel_replica_options(params).corr_time =
      ParametersLoadAccess::parallel_replica_options(params).corr_time_input /
      ParametersLoadAccess::constants(params).timeUnit;

  // [Temperature Accelerated Dynamics] //

  ParametersLoadAccess::tad_options(params).low_temperature =
      ini.GetReal("TAD", "low_temperature",
                  ParametersLoadAccess::tad_options(params).low_temperature);
  ParametersLoadAccess::tad_options(params).min_prefactor =
      ini.GetReal("TAD", "min_prefactor",
                  ParametersLoadAccess::tad_options(params).min_prefactor);
  ParametersLoadAccess::tad_options(params).confidence =
      ini.GetReal("TAD", "confidence",
                  ParametersLoadAccess::tad_options(params).confidence);

  // [Replica Exchange] //

  ParametersLoadAccess::replica_exchange_options(params)
      .temperature_distribution =
      toLowerCase(ini.Get("Replica Exchange", "temperature_distribution",
                          ParametersLoadAccess::replica_exchange_options(params)
                              .temperature_distribution));
  ParametersLoadAccess::replica_exchange_options(params).replicas =
      ini.GetInteger(
          "Replica Exchange", "replicas",
          ParametersLoadAccess::replica_exchange_options(params).replicas);
  ParametersLoadAccess::replica_exchange_options(params).exchange_trials =
      ini.GetInteger("Replica Exchange", "exchange_trials",
                     ParametersLoadAccess::replica_exchange_options(params)
                         .exchange_trials);
  ParametersLoadAccess::replica_exchange_options(params).sampling_time_input =
      ini.GetReal("Replica Exchange", "sampling_time",
                  ParametersLoadAccess::replica_exchange_options(params)
                      .sampling_time_input);
  ParametersLoadAccess::replica_exchange_options(params).sampling_time =
      ParametersLoadAccess::replica_exchange_options(params)
          .sampling_time_input /
      ParametersLoadAccess::constants(params).timeUnit;
  ParametersLoadAccess::replica_exchange_options(params).temperature_low =
      ini.GetReal("Replica Exchange", "temperature_low",
                  ParametersLoadAccess::main_options(params).temperature);
  ParametersLoadAccess::replica_exchange_options(params).temperature_high =
      ini.GetReal("Replica Exchange", "temperature_high",
                  ParametersLoadAccess::replica_exchange_options(params)
                      .temperature_high);
  ParametersLoadAccess::replica_exchange_options(params).exchange_period_input =
      ini.GetReal("Replica Exchange", "exchange_period",
                  ParametersLoadAccess::replica_exchange_options(params)
                      .exchange_period_input);
  ParametersLoadAccess::replica_exchange_options(params).exchange_period =
      ParametersLoadAccess::replica_exchange_options(params)
          .exchange_period_input /
      ParametersLoadAccess::constants(params).timeUnit;

  // [Hyperdynamics] //

  ParametersLoadAccess::hyperdynamics_options(params).rmd_time_input =
      ini.GetReal(
          "Hyperdynamics", "bb_rmd_time",
          ParametersLoadAccess::hyperdynamics_options(params).rmd_time_input);
  ParametersLoadAccess::hyperdynamics_options(params).rmd_time =
      ParametersLoadAccess::hyperdynamics_options(params).rmd_time_input /
      ParametersLoadAccess::constants(params).timeUnit;
  ParametersLoadAccess::hyperdynamics_options(params).boost_atom_list =
      toLowerCase(ini.Get(
          "Hyperdynamics", "bb_boost_atomlist",
          ParametersLoadAccess::hyperdynamics_options(params).boost_atom_list));
  ParametersLoadAccess::hyperdynamics_options(params).dvmax =
      ini.GetReal("Hyperdynamics", "bb_dvmax",
                  ParametersLoadAccess::hyperdynamics_options(params).dvmax);
  ParametersLoadAccess::hyperdynamics_options(params).qrr =
      ini.GetReal("Hyperdynamics", "bb_stretch_threshold",
                  ParametersLoadAccess::hyperdynamics_options(params).qrr);
  ParametersLoadAccess::hyperdynamics_options(params).prr =
      ini.GetReal("Hyperdynamics", "bb_ds_curvature",
                  ParametersLoadAccess::hyperdynamics_options(params).prr);
  ParametersLoadAccess::hyperdynamics_options(params).qcut =
      ini.GetReal("Hyperdynamics", "bb_rcut",
                  ParametersLoadAccess::hyperdynamics_options(params).qcut);
  ParametersLoadAccess::hyperdynamics_options(params).bias_potential =
      toLowerCase(ini.Get(
          "Hyperdynamics", "bias_potential",
          ParametersLoadAccess::hyperdynamics_options(params).bias_potential));

  // [Saddle Search] //

  ParametersLoadAccess::saddle_search_options(params).method = toLowerCase(
      ini.Get("Saddle Search", "method",
              ParametersLoadAccess::saddle_search_options(params).method));
  ParametersLoadAccess::saddle_search_options(params).minmode_method =
      toLowerCase(ini.Get(
          "Saddle Search", "min_mode_method",
          ParametersLoadAccess::saddle_search_options(params).minmode_method));
  ParametersLoadAccess::saddle_search_options(params).displace_magnitude =
      ini.GetReal("Saddle Search", "displace_magnitude",
                  ParametersLoadAccess::saddle_search_options(params)
                      .displace_magnitude);
  ParametersLoadAccess::saddle_search_options(params).displace_radius =
      ini.GetReal(
          "Saddle Search", "displace_radius",
          ParametersLoadAccess::saddle_search_options(params).displace_radius);
  ParametersLoadAccess::saddle_search_options(params).max_energy = ini.GetReal(
      "Saddle Search", "max_energy",
      ParametersLoadAccess::saddle_search_options(params).max_energy);
  ParametersLoadAccess::saddle_search_options(params).max_iterations =
      ini.GetInteger(
          "Saddle Search", "max_iterations",
          ParametersLoadAccess::optimizer_options(params).max_iterations);
  ParametersLoadAccess::saddle_search_options(params)
      .nonnegative_displacement_abort =
      ini.GetBoolean("Saddle Search", "nonnegative_displacement_abort",
                     ParametersLoadAccess::saddle_search_options(params)
                         .nonnegative_displacement_abort);
  ParametersLoadAccess::saddle_search_options(params).max_single_displace =
      ini.GetReal("Saddle Search", "max_single_displace",
                  ParametersLoadAccess::saddle_search_options(params)
                      .max_single_displace);
  ParametersLoadAccess::saddle_search_options(params).converged_force =
      ini.GetReal(
          "Saddle Search", "converged_force",
          ParametersLoadAccess::optimizer_options(params).converged_force);
  ParametersLoadAccess::saddle_search_options(params).perp_force_ratio =
      ini.GetReal(
          "Saddle Search", "perp_force_ratio",
          ParametersLoadAccess::saddle_search_options(params).perp_force_ratio);
  ParametersLoadAccess::saddle_search_options(params).displace_type =
      toLowerCase(ini.Get("Saddle Search", "client_displace_type",
                          eonc::EpiCenters::DISP_LOAD));
  ParametersLoadAccess::saddle_search_options(params).nonlocal_count_abort =
      ini.GetInteger("Saddle Search", "nonlocal_count_abort",
                     ParametersLoadAccess::saddle_search_options(params)
                         .nonlocal_count_abort);
  ParametersLoadAccess::saddle_search_options(params).nonlocal_distance_abort =
      ini.GetReal("Saddle Search", "nonlocal_distance_abort",
                  ParametersLoadAccess::saddle_search_options(params)
                      .nonlocal_distance_abort);
  if (params.saddle_search_options().displace_type !=
          eonc::EpiCenters::DISP_NOT_FCC_OR_HCP &&
      params.saddle_search_options().displace_type !=
          eonc::EpiCenters::DISP_MIN_COORDINATED &&
      params.saddle_search_options().displace_type !=
          eonc::EpiCenters::DISP_LAST_ATOM &&
      params.saddle_search_options().displace_type !=
          eonc::EpiCenters::DISP_RANDOM &&
      params.saddle_search_options().displace_type !=
          eonc::EpiCenters::DISP_LISTED_ATOMS) {
    ParametersLoadAccess::saddle_search_options(params).displace_type =
        eonc::EpiCenters::DISP_LOAD;
  }
  // Parse comma-separated atom list
  {
    std::string atomListStr =
        ini.Get("Saddle Search", "displace_atom_list", "");
    if (!atomListStr.empty()) {
      std::stringstream ss(atomListStr);
      std::string token;
      while (std::getline(ss, token, ',')) {
        size_t start = token.find_first_not_of(" \t");
        size_t end = token.find_last_not_of(" \t");
        if (start != std::string::npos) {
          ParametersLoadAccess::saddle_search_options(params)
              .displace_atom_list.push_back(
                  std::stol(token.substr(start, end - start + 1)));
        }
      }
    }
  }
  ParametersLoadAccess::saddle_search_options(params).confine_positive.enabled =
      ini.GetBoolean("Saddle Search", "confine_positive",
                     ParametersLoadAccess::saddle_search_options(params)
                         .confine_positive.enabled);
  if (ParametersLoadAccess::saddle_search_options(params)
          .confine_positive.enabled) {
    ParametersLoadAccess::saddle_search_options(params)
        .confine_positive.bowl_breakout =
        ini.GetBoolean("Saddle Search", "bowl_breakout",
                       ParametersLoadAccess::saddle_search_options(params)
                           .confine_positive.bowl_breakout);
    ParametersLoadAccess::saddle_search_options(params)
        .confine_positive.bowl_active =
        ini.GetInteger("Saddle Search", "bowl_active_atoms",
                       ParametersLoadAccess::saddle_search_options(params)
                           .confine_positive.bowl_active);
    ParametersLoadAccess::saddle_search_options(params)
        .confine_positive.min_force =
        ini.GetReal("Saddle Search", "confine_positive_min_move",
                    ParametersLoadAccess::saddle_search_options(params)
                        .confine_positive.min_force);
    ParametersLoadAccess::saddle_search_options(params)
        .confine_positive.scale_ratio =
        ini.GetReal("Saddle Search", "confine_positive_scale_ratio",
                    ParametersLoadAccess::saddle_search_options(params)
                        .confine_positive.scale_ratio);
    ParametersLoadAccess::saddle_search_options(params).confine_positive.boost =
        ini.GetReal("Saddle Search", "confine_positive_boost",
                    ParametersLoadAccess::saddle_search_options(params)
                        .confine_positive.boost);
    ParametersLoadAccess::saddle_search_options(params)
        .confine_positive.min_active =
        ini.GetInteger("Saddle Search", "confine_positive_min_active",
                       ParametersLoadAccess::saddle_search_options(params)
                           .confine_positive.min_active);
  }
  ParametersLoadAccess::saddle_search_options(params).dynamics.temperature =
      ParametersLoadAccess::main_options(params).temperature;
  ParametersLoadAccess::saddle_search_options(params).dynamics.temperature =
      ini.GetReal("Saddle Search", "dynamics_temperature",
                  ParametersLoadAccess::saddle_search_options(params)
                      .dynamics.temperature);
  ParametersLoadAccess::saddle_search_options(params)
      .dynamics.state_check_interval_input =
      ini.GetReal("Saddle Search", "dynamics_state_check_interval",
                  ParametersLoadAccess::saddle_search_options(params)
                      .dynamics.state_check_interval_input);
  ParametersLoadAccess::saddle_search_options(params)
      .dynamics.state_check_interval =
      ParametersLoadAccess::saddle_search_options(params)
          .dynamics.state_check_interval_input /
      ParametersLoadAccess::constants(params).timeUnit;
  ParametersLoadAccess::saddle_search_options(params)
      .dynamics.record_interval_input =
      ini.GetReal("Saddle Search", "dynamics_record_interval",
                  ParametersLoadAccess::saddle_search_options(params)
                      .dynamics.record_interval_input);
  ParametersLoadAccess::saddle_search_options(params).dynamics.record_interval =
      ParametersLoadAccess::saddle_search_options(params)
          .dynamics.record_interval_input /
      ParametersLoadAccess::constants(params).timeUnit;
  ParametersLoadAccess::saddle_search_options(params)
      .dynamics.linear_interpolation =
      ini.GetBoolean("Saddle Search", "dynamics_linear_interpolation",
                     ParametersLoadAccess::saddle_search_options(params)
                         .dynamics.linear_interpolation);
  ParametersLoadAccess::saddle_search_options(params).remove_rotation =
      ini.GetBoolean(
          "Saddle Search", "remove_rotation",
          ParametersLoadAccess::saddle_search_options(params).remove_rotation);
  ParametersLoadAccess::saddle_search_options(params)
      .dynamics.max_init_curvature =
      ini.GetReal("Saddle Search", "dynamics_max_init_curvature",
                  ParametersLoadAccess::saddle_search_options(params)
                      .dynamics.max_init_curvature);
  ParametersLoadAccess::saddle_search_options(params)
      .zero_mode_abort_curvature =
      ini.GetReal("Saddle Search", "zero_mode_abort_curvature",
                  ParametersLoadAccess::saddle_search_options(params)
                      .zero_mode_abort_curvature);

  // [Basin Hopping] //

  ParametersLoadAccess::basin_hopping_options(params).displacement =
      ini.GetReal(
          "Basin Hopping", "displacement",
          ParametersLoadAccess::basin_hopping_options(params).displacement);
  ParametersLoadAccess::basin_hopping_options(params).push_apart_distance =
      ini.GetReal("Basin Hopping", "push_apart_distance",
                  ParametersLoadAccess::basin_hopping_options(params)
                      .push_apart_distance);
  ParametersLoadAccess::basin_hopping_options(params)
      .initial_random_structure_probability =
      ini.GetReal("Basin Hopping", "initial_random_structure_probability",
                  ParametersLoadAccess::basin_hopping_options(params)
                      .initial_random_structure_probability);
  ParametersLoadAccess::basin_hopping_options(params).steps =
      ini.GetInteger("Basin Hopping", "steps",
                     ParametersLoadAccess::basin_hopping_options(params).steps);
  ParametersLoadAccess::basin_hopping_options(params).quenching_steps =
      ini.GetInteger(
          "Basin Hopping", "quenching_steps",
          ParametersLoadAccess::basin_hopping_options(params).quenching_steps);
  ParametersLoadAccess::basin_hopping_options(params).single_atom_displace =
      ini.GetBoolean("Basin Hopping", "single_atom_displace",
                     ParametersLoadAccess::basin_hopping_options(params)
                         .single_atom_displace);
  ParametersLoadAccess::basin_hopping_options(params).significant_structure =
      ini.GetBoolean("Basin Hopping", "significant_structure",
                     ParametersLoadAccess::basin_hopping_options(params)
                         .significant_structure);
  ParametersLoadAccess::basin_hopping_options(params).displacement_algorithm =
      toLowerCase(ini.Get("Basin Hopping", "displacement_algorithm",
                          ParametersLoadAccess::basin_hopping_options(params)
                              .displacement_algorithm));
  if (params.basin_hopping_options().displacement_algorithm != "standard" &&
      params.basin_hopping_options().displacement_algorithm != "linear" &&
      params.basin_hopping_options().displacement_algorithm != "quadratic") {
    EONC_LOG_ERROR("unknown displacement_algorithm {}",
                   ParametersLoadAccess::basin_hopping_options(params)
                       .displacement_algorithm);
    error = 1;
  }
  ParametersLoadAccess::basin_hopping_options(params)
      .displacement_distribution =
      toLowerCase(ini.Get("Basin Hopping", "displacement_distribution",
                          ParametersLoadAccess::basin_hopping_options(params)
                              .displacement_distribution));
  if (params.basin_hopping_options().displacement_distribution != "uniform" &&
      params.basin_hopping_options().displacement_distribution != "gaussian") {
    EONC_LOG_ERROR("unknown displacement_distribution {}",
                   ParametersLoadAccess::basin_hopping_options(params)
                       .displacement_distribution);
    error = 1;
  }
  ParametersLoadAccess::basin_hopping_options(params).swap_probability =
      ini.GetReal(
          "Basin Hopping", "swap_probability",
          ParametersLoadAccess::basin_hopping_options(params).swap_probability);
  ParametersLoadAccess::basin_hopping_options(params).jump_max = ini.GetInteger(
      "Basin Hopping", "jump_max",
      ParametersLoadAccess::basin_hopping_options(params).jump_max);
  ParametersLoadAccess::basin_hopping_options(params).jump_steps =
      ini.GetInteger(
          "Basin Hopping", "jump_steps",
          ParametersLoadAccess::basin_hopping_options(params).jump_steps);
  ParametersLoadAccess::basin_hopping_options(params).adjust_displacement =
      ini.GetBoolean("Basin Hopping", "adjust_displacement",
                     ParametersLoadAccess::basin_hopping_options(params)
                         .adjust_displacement);
  ParametersLoadAccess::basin_hopping_options(params).adjust_period =
      ini.GetInteger(
          "Basin Hopping", "adjust_period",
          ParametersLoadAccess::basin_hopping_options(params).adjust_period);
  ParametersLoadAccess::basin_hopping_options(params).adjust_fraction =
      ini.GetReal(
          "Basin Hopping", "adjust_fraction",
          ParametersLoadAccess::basin_hopping_options(params).adjust_fraction);
  ParametersLoadAccess::basin_hopping_options(params).target_ratio =
      ini.GetReal(
          "Basin Hopping", "target_ratio",
          ParametersLoadAccess::basin_hopping_options(params).target_ratio);
  ParametersLoadAccess::basin_hopping_options(params).write_unique =
      ini.GetBoolean(
          "Basin Hopping", "write_unique",
          ParametersLoadAccess::basin_hopping_options(params).write_unique);
  ParametersLoadAccess::basin_hopping_options(params).stop_energy = ini.GetReal(
      "Basin Hopping", "stop_energy",
      ParametersLoadAccess::basin_hopping_options(params).stop_energy);

  // [Global Optimization] //

  ParametersLoadAccess::global_optimization_options(params)
      .move_method = toLowerCase(ini.Get(
      "Global Optimization", "move_method",
      ParametersLoadAccess::global_optimization_options(params).move_method));
  ParametersLoadAccess::global_optimization_options(params).decision_method =
      toLowerCase(
          ini.Get("Global Optimization", "decision_method",
                  ParametersLoadAccess::global_optimization_options(params)
                      .decision_method));
  if (params.global_optimization_options().decision_method != "npew" &&
      params.global_optimization_options().decision_method != "boltzmann") {
    EONC_LOG_ERROR("unknown decision_method {}",
                   ParametersLoadAccess::global_optimization_options(params)
                       .decision_method);
    error = 1;
  }
  ParametersLoadAccess::global_optimization_options(params).steps =
      ini.GetInteger(
          "Global Optimization", "steps",
          ParametersLoadAccess::global_optimization_options(params).steps);
  ParametersLoadAccess::global_optimization_options(params).beta = ini.GetReal(
      "Global Optimization", "beta",
      ParametersLoadAccess::global_optimization_options(params).beta);
  ParametersLoadAccess::global_optimization_options(params).alpha = ini.GetReal(
      "Global Optimization", "alpha",
      ParametersLoadAccess::global_optimization_options(params).alpha);
  ParametersLoadAccess::global_optimization_options(params).mdmin =
      ini.GetInteger(
          "Global Optimization", "mdmin",
          ParametersLoadAccess::global_optimization_options(params).mdmin);
  ParametersLoadAccess::global_optimization_options(params).target_energy =
      ini.GetReal("Global Optimization", "target_energy",
                  ParametersLoadAccess::global_optimization_options(params)
                      .target_energy);

  // [BGSD] //

  ParametersLoadAccess::bgsd_options(params).alpha = ini.GetReal(
      "BGSD", "alpha", ParametersLoadAccess::bgsd_options(params).alpha);
  ParametersLoadAccess::bgsd_options(params).beta = ini.GetReal(
      "BGSD", "beta", ParametersLoadAccess::bgsd_options(params).beta);
  ParametersLoadAccess::bgsd_options(params).gradient_finite_difference =
      ini.GetReal("BGSD", "gradientfinitedifference",
                  ParametersLoadAccess::bgsd_options(params)
                      .gradient_finite_difference);
  ParametersLoadAccess::bgsd_options(params).grad2energy_convergence =
      ini.GetReal(
          "BGSD", "grad2energyconvergence",
          ParametersLoadAccess::bgsd_options(params).grad2energy_convergence);
  ParametersLoadAccess::bgsd_options(params).grad2force_convergence =
      ini.GetReal(
          "BGSD", "grad2forceconvergence",
          ParametersLoadAccess::bgsd_options(params).grad2force_convergence);

  // [Monte Carlo] //

  ParametersLoadAccess::monte_carlo_options(params).step_size =
      ini.GetReal("Monte Carlo", "step_size",
                  ParametersLoadAccess::monte_carlo_options(params).step_size);
  ParametersLoadAccess::monte_carlo_options(params).steps = static_cast<int>(
      ini.GetInteger("Monte Carlo", "steps",
                     ParametersLoadAccess::monte_carlo_options(params).steps));

  // [OH_TST] //

  ParametersLoadAccess::oh_tst_options(params).reactant_filename =
      ini.Get("OH_TST", "reactant_filename",
              ParametersLoadAccess::oh_tst_options(params).reactant_filename);
  ParametersLoadAccess::oh_tst_options(params).product_filename =
      ini.Get("OH_TST", "product_filename",
              ParametersLoadAccess::oh_tst_options(params).product_filename);
  ParametersLoadAccess::oh_tst_options(params).time_step =
      ini.GetReal("OH_TST", "time_step",
                  ParametersLoadAccess::oh_tst_options(params).time_step);
  ParametersLoadAccess::oh_tst_options(params).equil_steps =
      ini.GetInteger("OH_TST", "equil_steps",
                     ParametersLoadAccess::oh_tst_options(params).equil_steps);
  ParametersLoadAccess::oh_tst_options(params).sample_steps =
      ini.GetInteger("OH_TST", "sample_steps",
                     ParametersLoadAccess::oh_tst_options(params).sample_steps);
  ParametersLoadAccess::oh_tst_options(params).max_planes =
      ini.GetInteger("OH_TST", "max_planes",
                     ParametersLoadAccess::oh_tst_options(params).max_planes);
  ParametersLoadAccess::oh_tst_options(params).plane_mass =
      ini.GetReal("OH_TST", "plane_mass",
                  ParametersLoadAccess::oh_tst_options(params).plane_mass);
  ParametersLoadAccess::oh_tst_options(params).alpha_rot =
      ini.GetReal("OH_TST", "alpha_rot",
                  ParametersLoadAccess::oh_tst_options(params).alpha_rot);
  ParametersLoadAccess::oh_tst_options(params).plane_time_step =
      ini.GetReal("OH_TST", "plane_time_step",
                  ParametersLoadAccess::oh_tst_options(params).plane_time_step);
  ParametersLoadAccess::oh_tst_options(params).ds_max = ini.GetReal(
      "OH_TST", "ds_max", ParametersLoadAccess::oh_tst_options(params).ds_max);
  ParametersLoadAccess::oh_tst_options(params).dtheta_max =
      ini.GetReal("OH_TST", "dtheta_max",
                  ParametersLoadAccess::oh_tst_options(params).dtheta_max);
  ParametersLoadAccess::oh_tst_options(params).force_tol =
      ini.GetReal("OH_TST", "force_tol",
                  ParametersLoadAccess::oh_tst_options(params).force_tol);
  ParametersLoadAccess::oh_tst_options(params).s_init = ini.GetReal(
      "OH_TST", "s_init", ParametersLoadAccess::oh_tst_options(params).s_init);
  ParametersLoadAccess::oh_tst_options(params).reactant_md_steps =
      ini.GetInteger(
          "OH_TST", "reactant_md_steps",
          ParametersLoadAccess::oh_tst_options(params).reactant_md_steps);
  ParametersLoadAccess::oh_tst_options(params).symmetry_products =
      ini.Get("OH_TST", "symmetry_products",
              ParametersLoadAccess::oh_tst_options(params).symmetry_products);
  ParametersLoadAccess::oh_tst_options(params).max_delta_a =
      ini.GetReal("OH_TST", "max_delta_a",
                  ParametersLoadAccess::oh_tst_options(params).max_delta_a);
  ParametersLoadAccess::oh_tst_options(params).thermostat = toLowerCase(
      ini.Get("OH_TST", "thermostat",
              ParametersLoadAccess::oh_tst_options(params).thermostat));
  ParametersLoadAccess::oh_tst_options(params).gle_a_file =
      ini.Get("OH_TST", "gle_a_file",
              ParametersLoadAccess::oh_tst_options(params).gle_a_file);
  ParametersLoadAccess::oh_tst_options(params).pmf_scan =
      ini.GetBoolean("OH_TST", "pmf_scan",
                     ParametersLoadAccess::oh_tst_options(params).pmf_scan);
  ParametersLoadAccess::oh_tst_options(params).scan_planes =
      ini.GetInteger("OH_TST", "scan_planes",
                     ParametersLoadAccess::oh_tst_options(params).scan_planes);

  return error;
}

void validate_and_link(Parameters &params) {
  // Time unit conversions
  double tu = ParametersLoadAccess::constants(params).timeUnit;
  ParametersLoadAccess::optimizer_options(params).time_step =
      ParametersLoadAccess::optimizer_options(params).time_step_input / tu;
  ParametersLoadAccess::optimizer_options(params).max_time_step =
      ParametersLoadAccess::optimizer_options(params).max_time_step_input / tu;

  ParametersLoadAccess::dynamics_options(params).time_step =
      ParametersLoadAccess::dynamics_options(params).time_step_input / tu;
  ParametersLoadAccess::dynamics_options(params).time =
      ParametersLoadAccess::dynamics_options(params).time_input / tu;
  if (ParametersLoadAccess::dynamics_options(params).time_step > 0.0) {
    ParametersLoadAccess::dynamics_options(params).steps =
        static_cast<long>(std::floor(
            ParametersLoadAccess::dynamics_options(params).time /
                ParametersLoadAccess::dynamics_options(params).time_step +
            0.5));
  } else {
    ParametersLoadAccess::dynamics_options(params).steps = 0;
  }

  ParametersLoadAccess::thermostat_options(params).langevin_friction =
      ParametersLoadAccess::thermostat_options(params).langevin_friction_input *
      tu;

  ParametersLoadAccess::parallel_replica_options(params).dephase_time =
      ParametersLoadAccess::parallel_replica_options(params)
          .dephase_time_input /
      tu;
  ParametersLoadAccess::parallel_replica_options(params).state_check_interval =
      ParametersLoadAccess::parallel_replica_options(params)
          .state_check_interval_input /
      tu;
  ParametersLoadAccess::parallel_replica_options(params).record_interval =
      ParametersLoadAccess::parallel_replica_options(params)
          .record_interval_input /
      tu;
  ParametersLoadAccess::parallel_replica_options(params).corr_time =
      ParametersLoadAccess::parallel_replica_options(params).corr_time_input /
      tu;

  ParametersLoadAccess::replica_exchange_options(params).sampling_time =
      ParametersLoadAccess::replica_exchange_options(params)
          .sampling_time_input /
      tu;
  ParametersLoadAccess::replica_exchange_options(params).exchange_trials =
      ParametersLoadAccess::replica_exchange_options(params).replicas;

  ParametersLoadAccess::hyperdynamics_options(params).rmd_time =
      ParametersLoadAccess::hyperdynamics_options(params).rmd_time_input / tu;

  ParametersLoadAccess::saddle_search_options(params)
      .dynamics.state_check_interval =
      ParametersLoadAccess::saddle_search_options(params)
          .dynamics.state_check_interval_input /
      tu;
  ParametersLoadAccess::saddle_search_options(params).dynamics.record_interval =
      ParametersLoadAccess::saddle_search_options(params)
          .dynamics.record_interval_input /
      tu;

  // Cross-group defaults
  ParametersLoadAccess::process_search_options(params).minimization_offset =
      ParametersLoadAccess::optimizer_options(params).max_move;
  ParametersLoadAccess::neb_options(params).force_tolerance =
      ParametersLoadAccess::optimizer_options(params).converged_force;
}

} // namespace eonc::config
