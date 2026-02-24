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
#include "Matter.h"
#include "Parameters.h"
#include "Potential.h"
#include "version.h"

#ifdef WITH_SERVE_MODE
#include "ServeMode.h"
#endif

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>

using namespace std;

void singlePoint(std::unique_ptr<Matter> matter) {
  std::cout << "Energy:         " << std::fixed << std::setprecision(15)
            << matter->getPotentialEnergy() << std::endl;
  std::cout << "(free) Forces:         \n" << matter->getForcesFree() << "\n";
  std::cout << "Max atom force: " << std::scientific << matter->maxForce()
            << std::endl;
}

void minimize(std::unique_ptr<Matter> matter, const string &confileout) {
  matter->relax(false, false);
  if (!confileout.empty()) {
    std::cout << "Saving relaxed structure to " << confileout << std::endl;
  } else {
    std::cout << "No output file specified, not saving" << std::endl;
  }
  matter->matter2con(confileout);
}

void commandLine(int argc, char **argv) {
  bool sflag = false, mflag = false, pflag = false, cflag = false;
  double optConvergedForce = 0.001;
  string potential;
  string confile;
  string optimizer("cg");

  auto params = Parameters{};

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
      cxxopts::value<std::string>())
#ifdef WITH_SERVE_MODE
      ("serve",
       "Serve potential(s) over rgpot Cap'n Proto RPC. "
       "Spec: 'potential:port' or 'pot1:port1,pot2:port2'",
       cxxopts::value<std::string>())(
          "serve-host", "Host to bind RPC server(s) to",
          cxxopts::value<std::string>()->default_value("localhost"))(
          "serve-port", "Port for single-potential serve mode (used with -p)",
          cxxopts::value<uint16_t>()->default_value("12345"))(
          "replicas", "Number of replicated server instances (used with -p)",
          cxxopts::value<size_t>()->default_value("1"))(
          "gateway",
          "Run a single gateway port backed by N pool instances "
          "(use with -p and --replicas)",
          cxxopts::value<bool>()->default_value("false"))(
          "config",
          "Config file for potential parameters (INI format, "
          "e.g. [Metatomic] model_path=model.pt)",
          cxxopts::value<std::string>())
#endif
          ("h,help", "Print usage");

  try {
    auto result = options.parse(argc, argv);

    if (result.count("help")) {
      std::cout << options.help() << std::endl;
      exit(0);
    }

    if (result.count("version")) {
      std::cout << "eonclient version " << VERSION << " (" << GIT_HASH << ")"
                << std::endl;
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

    if (result.count("tolerance")) {
      params.structure_comparison_options.distance_difference =
          result["tolerance"].as<double>();
    }

    if (sflag && mflag) {
      std::cerr << "Cannot specify both minimization and single point"
                << std::endl;
      exit(2);
    }

    if (!pflag && (sflag || mflag)) {
      std::cerr << "Must specify a potential" << std::endl;
      exit(2);
    }

#ifdef WITH_SERVE_MODE
    // Load config file if provided (for potential-specific parameters
    // like model_path, device, length_unit, etc.)
    if (result.count("config")) {
      auto config_path = result["config"].as<std::string>();
      std::ifstream config_file(config_path);
      if (!config_file.is_open()) {
        std::cerr << "Cannot open config file: " << config_path << std::endl;
        exit(2);
      }
      params.load(config_path);
    }

    // Handle --serve mode (does not require a con file)
    if (result.count("serve")) {
      auto spec = result["serve"].as<std::string>();
      auto endpoints = parseServeSpec(spec);
      if (endpoints.empty()) {
        std::cerr << "No valid serve endpoints in spec: " << spec << std::endl;
        exit(2);
      }
      serveMultiple(endpoints, params);
      exit(0);
    }

    // Handle -p with serve flags (single potential serve mode)
    if (pflag && !sflag && !mflag && !cflag &&
        (result.count("serve-port") || result.count("replicas") ||
         result.count("gateway"))) {
      for (auto &ch : potential) {
        ch = tolower(ch);
      }
      params.potential_options.potential =
          magic_enum::enum_cast<PotType>(potential,
                                         magic_enum::case_insensitive)
              .value_or(PotType::UNKNOWN);
      auto host = result["serve-host"].as<std::string>();
      auto port = result["serve-port"].as<uint16_t>();
      auto reps = result["replicas"].as<size_t>();
      bool gw = result["gateway"].as<bool>();

      if (gw) {
        serveGateway(params, host, port, reps);
      } else if (reps > 1) {
        serveReplicated(params, host, port, reps);
      } else {
        serveMode(params, host, port);
      }
      exit(0);
    }

    // Config-driven serve (no -p or --serve, just --config with [Serve])
    if (!pflag && !sflag && !mflag && !cflag && result.count("config") &&
        !result.count("serve") &&
        (!params.serve_options.endpoints.empty() ||
         params.serve_options.gateway_port > 0 ||
         params.serve_options.replicas > 1)) {
      serveFromConfig(params);
      exit(0);
    }
#endif

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

    if (!cflag) {
      params.potential_options.potential =
          magic_enum::enum_cast<PotType>(potential,
                                         magic_enum::case_insensitive)
              .value_or(PotType::UNKNOWN);
    }

    if (!sflag) {
      params.optimizer_options.method =
          magic_enum::enum_cast<OptType>(optimizer,
                                         magic_enum::case_insensitive)
              .value_or(OptType::CG);
      params.optimizer_options.converged_force = optConvergedForce;
    }

    auto pot = helper_functions::makePotential(params);
    auto matter = std::make_unique<Matter>(pot, params);
    auto matter2 = std::make_unique<Matter>(pot, params);
    matter->con2matter(confile);

    string confileout;
    if (unmatched.size() == 2) {
      confileout = unmatched[1];
      if (cflag)
        matter2->con2matter(confileout);
    }

    if (sflag) {
      singlePoint(std::move(matter));
    } else if (mflag) {
      minimize(std::move(matter), confileout);
    } else if (cflag) {
      params.structure_comparison_options.check_rotation = true;
      if (matter->compare(*matter2, true)) {
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
