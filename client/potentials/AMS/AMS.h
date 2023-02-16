//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef AMS_POT
#define AMS_POT

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
#include <fmt/core.h>
#include <fmt/format.h>
#include <fstream>
#include <string>
#include <filesystem>

class AMS : public Potential {

public:
  AMS(Parameters *p);
  ~AMS();
  void initialize(){};
  void cleanMemory(void);
  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, const double *box);
  std::pair<double, AtomMatrix> get_ef(const AtomMatrix pos,
                                       const VectorXi atmnrs,
                                       const Matrix3d box) override {
    double energy{std::numeric_limits<double>::infinity()};
    long nAtoms{pos.rows()};
    AtomMatrix forces{Eigen::MatrixXd::Zero(nAtoms, 3)};
    this->force(nAtoms, pos.data(), atmnrs.data(), forces.data(), &energy,
                box.data());
    return std::make_pair(energy, forces);
  };

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
  std::string generate_run(Parameters *p);
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

#endif
