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
#include "CommandLine.h"
#include "Job.h"
#include "Parameters.h"
#include "Potential.h"
#include "client/io/WriteCreator.hpp"
#include "client/matter/Matter.h"
#include "client/matter/MatterCreator.hpp"
#include "client/matter/StructComparer.hpp"
#include "version.h"

#include <cstdlib>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
namespace eonc {
using namespace std;

void minimize(std::unique_ptr<Matter> matter, const string &confileout) {
  // XXX: Fix this
  // matter->relax(false, false);
  if (!confileout.empty()) {
    std::cout << "Saving relaxed structure to " << confileout << std::endl;
  } else {
    std::cout << "No output file specified, not saving" << std::endl;
  }
  const auto config = toml::table{{"Main", toml::table{{"write", "con"}}}};
  auto con = eonc::io::mkWriter(config);
  con->write(*matter, confileout);
}

void commandLine(int argc, char **argv) {
  bool sflag = false, mflag = false, pflag = false, cflag = false;
  double optConvergedForce = 0.001;
  string potential;
  string confile;
  string optimizer("cg");

  auto params = std::make_shared<Parameters>();

  cxxopts::Options options("eonclient", "The eOn client");
  options.add_options()("v,version", "Print version information")(
      "m,minimize", "Minimization of inputConfile saves to outputConfile")(
      "s,single", "Single point energy of inputConfile")(
      "c,compare", "Compare structures of inputConfile to outputConfile")(
      "o,optimizer", "Optimization method",
      cxxopts::value<std::string>()->default_value("cg"))(
      "f,force", "Convergence force",
      cxxopts::value<double>()->default_value("0.001"))(
      "t,tolerance", "Distance tolerance",
      cxxopts::value<double>()->default_value("0.1"))(
      "p,potential", "The potential (e.g. qsc, lj, eam_al)",
      cxxopts::value<std::string>())("h,help", "Print usage");

  try {
    auto result = options.parse(argc, argv);

    if (result.count("help")) {
      std::cout << options.help() << std::endl;
      exit(0);
    }

    if (result.count("version")) {
      std::cout << "eonclient version r" << VERSION << std::endl;
      std::cout << "          compiled " << BUILD_DATE << std::endl;
      exit(0);
    }

    if (result.count("minimize")) {
      mflag = true;
    }

    if (result.count("single")) {
      sflag = true;
    }

    if (result.count("compare")) {
      cflag = true;
    }

    if (result.count("potential")) {
      pflag = true;
      potential = result["potential"].as<std::string>();
    }

    if (result.count("optimizer")) {
      optimizer = result["optimizer"].as<std::string>();
    }

    if (result.count("force")) {
      optConvergedForce = result["force"].as<double>();
    }

    // TODO(rg):: Use this
    // if (result.count("tolerance")) {
    //   params->structcomp.distanceDifference =
    //   result["tolerance"].as<double>();
    // }

    if (sflag && mflag) {
      std::cerr << "Cannot specify both minimization and single point"
                << std::endl;
      exit(2);
    }

    if (!pflag && (sflag || mflag)) {
      std::cerr << "Must specify a potential" << std::endl;
      exit(2);
    }

    if (!cflag) {
      for (auto &ch : potential) {
        ch = tolower(ch);
      }
    }

    auto unmatched = result.unmatched();
    if (unmatched.size() < 1) {
      std::cerr << "At least one non-option argument is required: the con file"
                << std::endl;
      exit(2);
    } else {
      confile = unmatched[0];
    }

    auto tbl =
        toml::table{{"Potential", toml::table{{"potential", "unknown"}}}};
    if (!cflag) {
      tbl["Potential"].as_table()->insert_or_assign("potential", potential);
    } else {
      // For structure comparison, we need only to ensure a potential exists..
      tbl["Potential"].as_table()->insert_or_assign("potential", "lj");
    }

    if (!sflag) {
      params->optim.method = magic_enum::enum_cast<OptType>(
                                 optimizer, magic_enum::case_insensitive)
                                 .value_or(OptType::CG);
      params->optim.convergedForce = optConvergedForce;
    }

    auto pot = makePotential(tbl);
    auto mat1 = Matter(pot);
    eonc::mat::from_con(mat1, confile);

    string confileout;
    if (unmatched.size() == 2) {
      confileout = unmatched[1];
    }

    if (sflag) {
      auto tbl = toml::table{{"Main", toml::table{{"job", "point"}}}};
      // Run PointJob
      auto spj = eonc::makeJob(tbl, mat1);
      auto res = spj->run();
    } else if (mflag) {
      // XXX: Finish
      // minimize(matter, confileout);
    } else if (cflag) {
      if (unmatched.size() != 2) {
        throw std::runtime_error("Comparison needs two files!");
      }
      auto mat2 = Matter(pot);
      eonc::mat::from_con(mat2, confileout);
      auto tbl = toml::table{
          {"Structure_Comparison", toml::table{{"checkRotation", true}}}};
      auto sc = eonc::mat::StructComparer(tbl);
      if (sc.compare(mat1, mat2, true)) {
        std::cout << "Structures match" << std::endl;
      } else {
        std::cout << "Structures do not match" << std::endl;
      }
    }
  } catch (const cxxopts::exceptions::exception &e) {
    std::cerr << "Error parsing options: " << e.what() << std::endl;
    std::cerr << options.help() << std::endl;
    exit(1);
  }
}
} // namespace eonc
