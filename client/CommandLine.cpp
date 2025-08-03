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
#include "client/Parser.hpp"
#include "version.h"

#include <cstdlib>
#include <iostream>
#include <memory>
#include <spdlog/spdlog.h>
#include <string>

#include "thirdparty/cxxopts.hpp"
namespace eonc {
using namespace std;

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
      "input", "Input files",
      cxxopts::value<std::vector<std::string>>()->implicit_value({"pos.con"}))(
      "o,output", "Output file", cxxopts::value<std::string>())(
      "opt", "Optimization method",
      cxxopts::value<std::string>()->default_value("cg"))(
      "f,force", "Convergence force",
      cxxopts::value<double>()->default_value("0.001"))(
      "metric", "Convergence metric",
      cxxopts::value<std::string>()->default_value("norm"))(
      "t,tolerance", "Distance tolerance",
      cxxopts::value<double>()->default_value("0.1"))(
      "check_rotation", "Check rotation for structure comparison",
      cxxopts::value<bool>()->default_value("true"))(
      "p,potential", "The potential (e.g. qsc, lj, eam_al)",
      cxxopts::value<std::string>())(
      "q,quiet", "Pure CLI mode, no version and timing info",
      cxxopts::value<bool>()->implicit_value("false"))(
      "conf", "Configuration file",
      cxxopts::value<std::string>()->default_value("config.toml"))(
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

    // SPDLOG_LOGGER_DEBUG(log, argc);
    if (result.count("conf") || (argc == 1)) {
      std::cout << "Loading " << result["conf"].as<std::string>()
                << " and overwriting with commandline options" << std::endl;
      tbl = eonc::loadTOML(result["conf"].as<std::string>());

      // Check if 'inputs' is present in the TOML file; if not, assign a default
      if (!tbl["Main"].as_table()->contains("inputs")) {
        SPDLOG_LOGGER_WARN(log, "Input file key missing, assuming pos.con");
        tbl["Main"].as_table()->insert_or_assign("inputs",
                                                 toml::array{"pos.con"});
      }
    }

    if (result.count("minimize")) {
      tbl["Main"].as_table()->insert_or_assign("job", "minimization");
      // Constructs a default optimizer setup
      tbl.insert_or_assign("Optimizer", toml::table{});
      tbl["Optimizer"].as_table()->insert_or_assign(
          "method", result["opt"].as<std::string>());
      tbl["Optimizer"].as_table()->insert_or_assign(
          "converged_force", result["force"].as<double>());
      tbl["Optimizer"].as_table()->insert_or_assign(
          "convergence_metric", result["metric"].as<std::string>());
    }

    if (result.count("compare")) {
      tbl["Main"].as_table()->insert_or_assign("job", "structure_comparison");
    }

    if (result.count("potential")) {
      tbl["Potential"].as_table()->insert_or_assign(
          "potential", result["potential"].as<std::string>());
    }

    if (result.count("opt")) {
      tbl.insert_or_assign("Optimizer", toml::table{});
      tbl["Optimizer"].as_table()->insert_or_assign(
          "method", result["opt"].as<std::string>());
      tbl["Optimizer"].as_table()->insert_or_assign(
          "converged_force", result["force"].as<double>());
      tbl["Optimizer"].as_table()->insert_or_assign(
          "convergence_metric", result["metric"].as<std::string>());
    }

    if (result.count("input")) {
      toml::array inputs_array;
      for (const auto &input : result["input"].as<std::vector<std::string>>()) {
        inputs_array.push_back(input);
      }
      tbl["Main"].as_table()->insert_or_assign("inputs",
                                               std::move(inputs_array));
    }

    // TODO(rg):: Use this
    // if (result.count("tolerance")) {
    //   params->structcomp.distanceDifference =
    //   result["tolerance"].as<double>();
    // }

    if (result.count("quiet")) {
      tbl["Main"].as_table()->insert_or_assign("quiet",
                                               result["quiet"].as<bool>());
    }

    if (result.count("output")) {
      tbl["Main"].as_table()->insert_or_assign(
          "output", result["output"].as<std::string>());
    }

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

    if (result.count("compare")) {
      // For structure comparison, we need only to ensure a potential exists..
      tbl["Potential"].as_table()->insert_or_assign("potential", "lj");
      tbl.insert_or_assign("Structure_Comparison", toml::table{});
      if (result.count("check_rotation")) {
        tbl["Structure_Comparison"].as_table()->insert_or_assign(
            "check_rotation", result["check_rotation"].as<bool>());
      }
      if (result.count("tolerance")) {
        tbl["Structure_Comparison"].as_table()->insert_or_assign(
            "distance_difference", result["tolerance"].as<double>());
      }
    }

    // Assign defaults if not present so far..
    // Check if 'output' is present in the TOML file; if not, assign a default
    if (!tbl["Main"].as_table()->contains("output")) {
      SPDLOG_LOGGER_WARN(log, "Output key missing, assuming res.con");
      tbl["Main"].as_table()->insert_or_assign("output",
                                               toml::array{"res.con"});
    }

  } catch (const cxxopts::exceptions::exception &e) {
    std::cerr << "Error parsing options: " << e.what() << std::endl;
    std::cerr << options.help() << std::endl;
    exit(EXIT_FAILURE);
  }
  std::cout << tbl << std::endl << std::endl;
  return tbl;
}
} // namespace eonc
