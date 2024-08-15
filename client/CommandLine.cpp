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
#include "Job.hpp"
#include "client/Parser.hpp"
#include "client/io/WriteCreator.hpp"
#include "client/matter/Matter.h"
#include "version.h"

#include <cstdlib>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

#include "thirdparty/cxxopts.hpp"
namespace eonc {
using namespace std;

void minimize(Matter &matter, const string &confileout) {
  auto tbl = toml::table{
      {"Main", toml::table{{"job", "minimization"}, {"write", "con"}}}};
  auto relax = mkJob(tbl);
  bool result = JobRunner(relax, matter);
  if (!result) {
    throw std::runtime_error("Something went wrong running the job");
  }
  if (!confileout.empty()) {
    std::cout << "Saving relaxed structure to " << confileout << std::endl;
  } else {
    std::cout << "No output file specified, not saving" << std::endl;
  }
  auto con = eonc::io::mkWriter(tbl);
  con->write(matter, confileout);
}

toml::table commandLine(std::shared_ptr<spdlog::logger> log, int argc,
                        char **argv) {
  // Default configuration
  auto tbl = toml::table{{
                             "Main",
                             toml::table{{"job", "point"}, {"write", "con"}},
                         },
                         {"Potential", toml::table{{"potential", "unknown"}}}};

  cxxopts::Options options("eonclient", "The eOn client");
  options.add_options()("v,version", "Print version information")(
      "m,minimize", "Minimization of inputConfile saves to outputConfile")(
      "s,single", "Single point energy of inputConfile")(
      "c,compare", "Compare structures of inputConfile to outputConfile")(
      "i,input", "Input files",
      cxxopts::value<std::vector<std::string>>()->implicit_value({"pos.con"}))(
      "o,output", "Output file", cxxopts::value<std::string>())(
      "opt", "Optimization method",
      cxxopts::value<std::string>()->implicit_value("cg"))(
      "f,force", "Convergence force",
      cxxopts::value<double>()->default_value("0.001"))(
      "t,tolerance", "Distance tolerance",
      cxxopts::value<double>()->default_value("0.1"))(
      "p,potential", "The potential (e.g. qsc, lj, eam_al)",
      cxxopts::value<std::string>())(
      "conf", "Configuration file",
      cxxopts::value<std::string>()->implicit_value("config.toml"))(
      "h,help", "Print usage");

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

    if (result.count("config")) {
      std::cout << "Loading " << result["conf"].as<std::string>()
                << " and overwriting with commandline options";
      tbl = eonc::loadTOML(result["conf"].as<std::string>());
    } else {
      // For suppresssing version and timing information
      tbl["Main"].as_table()->insert_or_assign("suppress", true);
    }

    if (result.count("minimize")) {
      tbl["Main"].as_table()->insert_or_assign("job", "minimize");
    }

    if (result.count("compare")) {
      tbl["Main"].as_table()->insert_or_assign("job", "structure_compare");
    }

    if (result.count("potential")) {
      tbl["Potential"].as_table()->insert_or_assign(
          "potential", result["potential"].as<std::string>());
    }

    if (result.count("opt")) {
      tbl.insert_or_assign("Optimizer", toml::table{});
      tbl["Optimizer"].as_table()->insert_or_assign(
          "opt_method", result["opt"].as<std::string>());
      tbl["Optimizer"].as_table()->insert_or_assign(
          "converged_force", result["force"].as<double>());
    }

    if (result.count("input")) {
      toml::array inputs_array;
      for (const auto &input : result["input"].as<std::vector<std::string>>()) {
        inputs_array.push_back(input);
      }
      tbl["Main"].as_table()->insert_or_assign("inputs",
                                               std::move(inputs_array));
    }

    if (result.count("output")) {
      tbl["Main"].as_table()->insert_or_assign(
          "output", result["output"].as<std::string>());
    }

    // TODO(rg):: Use this
    // if (result.count("tolerance")) {
    //   params->structcomp.distanceDifference =
    //   result["tolerance"].as<double>();
    // }

    if (result.count("point") && result.count("minimize")) {
      std::cerr << "Cannot specify both minimization and single point"
                << std::endl;
      exit(EXIT_FAILURE);
    }

    if (!result.count("potential") &&
        (result.count("point") || result.count("minimize"))) {
      std::cerr << "Must specify a potential" << std::endl;
      exit(EXIT_FAILURE);
    }

    // auto unmatched = result.unmatched();
    // if (unmatched.size() < 1) {
    //   std::cerr << "At least one non-option argument is required: the con
    //   file"
    //             << std::endl;
    //   exit(EXIT_FAILURE);
    // } else {
    //   tbl["Main"].as_table()->insert_or_assign("ep_one", unmatched[0]);
    // }

    // std::string confile = tbl["Main"]["ep_one"].value<std::string>().value();

    if (result.count("compare")) {
      // For structure comparison, we need only to ensure a potential exists..
      tbl["Potential"].as_table()->insert_or_assign("potential", "lj");
    }

    // TODO(rg) Port
    // if (!sflag) {
    //   params->optim.method = magic_enum::enum_cast<OptType>(
    //                              optimizer, magic_enum::case_insensitive)
    //                              .value_or(OptType::CG);
    //   params->optim.convergedForce = optConvergedForce;
    // }

    // auto CACHELOT_CMD = cachelot::cache::Cache::Create(
    //     eonc::cache::cache_memory, eonc::cache::page_size,
    //     eonc::cache::hash_initial, true);
    // auto pcache = eonc::cache::PotentialCache();
    // pcache.set_cache(&CACHELOT_CMD);

    // auto pot = makePotential(tbl);
    // pot->set_cache(&pcache);

    // auto mat1 = Matter(pot);
    // eonc::mat::ConFileParser cfp;
    // cfp.parse(mat1, confile);

    // string confileout;
    // if (unmatched.size() == 2) {
    //   tbl["Main"].as_table()->insert_or_assign("ep_two", unmatched[1]);
    // }

    // if (sflag) {
    //   auto tbl = toml::table{{"Main", toml::table{{"job", "point"}}}};
    //   // Run PointJob
    //   auto spj = mkJob(tbl);
    //   bool result = JobRunner(spj, mat1);
    //   if (!result) {
    //     throw std::runtime_error("Something went wrong running the job");
    //   }
    // } else if (mflag) {
    //   minimize(mat1, confileout);
    // } else if (cflag) {
    //   if (unmatched.size() != 2) {
    //     throw std::runtime_error("Comparison needs two files!");
    //   }
    //   auto mat2 = Matter(pot);
    //   cfp.parse(mat2, confileout);
    //   auto tbl = toml::table{
    //       {"Structure_Comparison", toml::table{{"checkRotation", true}}}};
    //   auto sc = eonc::mat::StructComparer(tbl);
    //   if (sc.compare(mat1, mat2, true)) {
    //     std::cout << "Structures match" << std::endl;
    //   } else {
    //     std::cout << "Structures do not match" << std::endl;
    //   }
    // }
  } catch (const cxxopts::exceptions::exception &e) {
    std::cerr << "Error parsing options: " << e.what() << std::endl;
    std::cerr << options.help() << std::endl;
    exit(EXIT_FAILURE);
  }
  // std::cout << tbl;
  return tbl;
}
} // namespace eonc
