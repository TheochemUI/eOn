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

#pragma once
#include "../../Matter.h"
#include "../../Potential.h"

#include <absl/strings/numbers.h>
#include <absl/strings/str_cat.h>
#include <absl/strings/str_join.h>
#include <absl/strings/str_split.h>
#include <absl/strings/string_view.h>
#include <boost/asio.hpp>
#include <boost/process.hpp>
#include <boost/process/environment.hpp>

#include <algorithm>
#include <filesystem>
#include <fmt/core.h>
#include <fmt/format.h>
#include <fstream>
#include <string>

class AMS : public Potential {

public:
  AMS(Parameters& p);
  ~AMS();
  void initialize(){};
  void cleanMemory(void);
  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box);

private:
  //!< Creates a script to run AMS
  void passToSystem(long N, const double *R, const int *atomicNrs,
                    const double *box);
  void smallSys(long N, const double *R, const int *atomicNrs,
                const double *box);
  // Switch between jobs
  void switchjob();
  // Write restart files
  void write_restart();
  std::string engine, forcefield, model, xc, resources, basis;
  std::string engine_setup, engine_lower;
  // Generate run configuration
  std::string generate_run(Parameters& p);
  // Environment
  boost::process::native_environment nativenv;
  int amsevals;
  bool first_run, can_restart;
  std::string cjob, pjob;
  std::ofstream restartFrom;
  const double forceConversion =
      -51.4220862; // Forces from hartree/bohr to eV/Angstrom, -1 for the
                   // gradients
  const double energyConversion = 27.2114; // Energy in hartree to eV
  const double lengthConversion = 1.88973; // Coordinates from Angstrom to Bohr
  void runAMS();
  void updateCoord(long N, const double *R);
  void extract_rkf(long N, std::string key);
  double extract_scalar_rkf(std::string key); // single quantity
  std::vector<double>
  extract_cartesian_rkf(std::string key); // 3 x N quantities
  // Debugging utilities
  bool validate_order();
  std::string readFile(std::filesystem::path path);
  void recieveFromSystem(long N, double *F, double *U);
};
