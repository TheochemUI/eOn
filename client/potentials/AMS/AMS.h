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

#include <algorithm>
#include <fstream>
#include <string>

class AMS : public Potential {

public:
  AMS(Parameters *p);
  ~AMS();
  void initialize(){};
  void cleanMemory(void);
  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, const double *box, int nImages);

private:
  //!< Creates a script to run AMS
  void passToSystem(long N, const double *R, const int *atomicNrs,
                    const double *box);
  void recieveFromSystem(long N, double *F, double *U);
  const char *engine;
  const char *model;
  const char *forcefield;
  const char *xc;
  int counter;
  bool first_run, job_one;
  boost::asio::io_context ams_run, ams_rkf;
  std::future<std::string> err, errg, erre, errc, errec, ecdump, edump, gdump,
      cdump;
  std::string jname, restartj, coordDump;
  std::ofstream updCoord, restartFrom;
  double x, energy;
  std::vector<double> gradients;
  std::vector<std::string> tmp, energ, grad;
  const double forceConversion =
      -51.4220862; // Forces from hartree/bohr to eV/Angstrom
  const double energyConversion = 27.2114; // Energy in hartree to eV
};

#endif
