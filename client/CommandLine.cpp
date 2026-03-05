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

#include <argum.h>

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>

using namespace Argum;
using namespace std;

// Create a colorizer with default color scheme
constexpr auto colorScheme = basicDefaultColorScheme<char>;
BasicColorizer<char> colorizer(colorScheme);

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

void printFeatures() {
  std::cout << colorizer.heading("Compile-time features:") << std::endl;
  // Parse FEATURES_STRING and colorize enabled/disabled
  std::istringstream stream(FEATURES_STRING);
  std::string line;
  while (std::getline(stream, line)) {
    if (line.find(": enabled") != std::string::npos) {
      // Extract feature name and colorize "enabled" with warning (yellow)
      size_t colonPos = line.find(": ");
      if (colonPos != std::string::npos) {
        std::string featureName = line.substr(0, colonPos + 2);
        std::string status = line.substr(colonPos + 2);
        std::cout << featureName << colorizer.warning(status) << std::endl;
      } else {
        std::cout << line << std::endl;
      }
    } else if (line.find(": disabled") != std::string::npos) {
      // Extract feature name and colorize "disabled" with error (red)
      size_t colonPos = line.find(": ");
      if (colonPos != std::string::npos) {
        std::string featureName = line.substr(0, colonPos + 2);
        std::string status = line.substr(colonPos + 2);
        std::cout << featureName << colorizer.error(status) << std::endl;
      } else {
        std::cout << line << std::endl;
      }
    } else {
      std::cout << line << std::endl;
    }
  }
}

void commandLine(int argc, char **argv) {
  bool sflag = false, mflag = false, pflag = false, cflag = false;
  double optConvergedForce = 0.001;
  string potential;
  string confile;
  string optimizer("cg");
  optional<string> config_path;

#ifdef WITH_SERVE_MODE
  optional<string> serve_spec;
  optional<string> serve_host("localhost");
  optional<uint16_t> serve_port(12345);
  optional<size_t> replicas(1);
  optional<bool> gateway(false);
#endif

  auto params = Parameters{};

  const char *progname = (argc ? argv[0] : "eonclient");

  Parser parser;

  parser.add(Option("--help", "-h")
                 .help("show this help message and exit")
                 .handler([&]() {
                   // Format help with color
                   auto helpText = parser.formatHelp(progname);
                   cout << colorizer.heading("eOn Client - Help") << "\n\n";
                   cout << helpText;
                   exit(EXIT_SUCCESS);
                 }));

  parser.add(Option("--version", "-v")
                 .help("Print version information")
                 .handler([&]() {
                   cout << VERSION_STRING << endl;
                   exit(EXIT_SUCCESS);
                 }));

  parser.add(
      Option("--features").help("Print compile-time features").handler([&]() {
        printFeatures();
        exit(EXIT_SUCCESS);
      }));

  parser.add(Option("--minimize", "-m")
                 .help("Minimization of inputConfile saves to outputConfile")
                 .handler([&]() { mflag = true; }));

  parser.add(Option("--single", "-s")
                 .help("Single point energy of inputConfile")
                 .handler([&]() { sflag = true; }));

  parser.add(Option("--compare", "-c")
                 .help("Compare structures of inputConfile to outputConfile")
                 .handler([&]() { cflag = true; }));

  parser.add(
      Option("--optimizer", "-o")
          .argName("METHOD")
          .help("Optimization method")
          .handler([&](const string_view &value) { optimizer = value; }));

  parser.add(Option("--force", "-f")
                 .argName("VALUE")
                 .help("Convergence force")
                 .handler([&](const string_view &value) {
                   optConvergedForce = parseFloatingPoint<double>(value);
                 }));

  parser.add(Option("--tolerance", "-t")
                 .argName("VALUE")
                 .help("Distance tolerance")
                 .handler([&](const string_view &value) {
                   params.structure_comparison_options.distance_difference =
                       parseFloatingPoint<double>(value);
                 }));

  parser.add(Option("--potential", "-p")
                 .argName("POTENTIAL")
                 .help("The potential (e.g. qsc, lj, eam_al)")
                 .handler([&](const string_view &value) {
                   pflag = true;
                   potential = value;
                 }));

#ifdef WITH_SERVE_MODE
  parser.add(
      Option("--serve")
          .argName("SPEC")
          .help("Serve potential(s) over rgpot Cap'n Proto RPC. "
                "Spec: 'potential:port' or 'pot1:port1,pot2:port2'")
          .handler([&](const string_view &value) { serve_spec = value; }));

  parser.add(
      Option("--serve-host")
          .argName("HOST")
          .help("Host to bind RPC server(s) to")
          .handler([&](const string_view &value) { serve_host = value; }));

  parser.add(Option("--serve-port")
                 .argName("PORT")
                 .help("Port for single-potential serve mode (used with -p)")
                 .handler([&](const string_view &value) {
                   serve_port = parseIntegral<uint16_t>(value);
                 }));

  parser.add(Option("--replicas")
                 .argName("N")
                 .help("Number of replicated server instances (used with -p)")
                 .handler([&](const string_view &value) {
                   replicas = parseIntegral<size_t>(value);
                 }));

  parser.add(Option("--gateway")
                 .help("Run a single gateway port backed by N pool instances "
                       "(use with -p and --replicas)")
                 .handler([&]() { gateway = true; }));

  parser.add(
      Option("--config")
          .argName("FILE")
          .help("Config file for potential parameters (INI format, "
                "e.g. [Metatomic] model_path=model.pt)")
          .handler([&](const string_view &value) { config_path = value; }));
#endif

  parser.add(Positional("confile")
                 .help("Input structure file")
                 .occurs(zeroOrMoreTimes)
                 .handler([&](const string_view &value) { confile = value; }));

  parser.add(Positional("confileout")
                 .help("Output structure file (optional)")
                 .occurs(zeroOrMoreTimes)
                 .handler([&](const string_view &value) {
                   // confileout will be set after confile
                 }));

  try {
    parser.parse(argc, argv);
  } catch (const ParsingException &ex) {
    cerr << colorizer.error(ex.message()) << '\n';
    cerr << colorizer.warning(parser.formatUsage(progname)) << '\n';
    exit(EXIT_FAILURE);
  }

  // Re-parse positional arguments since argum handles them in order
  // We need to manually extract confile and confileout
  // Reset and re-parse for positional args
  confile.clear();
  string confileout;
  int positional_count = 0;
  for (int i = 1; i < argc; ++i) {
    string arg = argv[i];
    if (arg.starts_with("-")) {
      continue;
    }
    if (positional_count == 0) {
      confile = arg;
    } else if (positional_count == 1) {
      confileout = arg;
    }
    ++positional_count;
  }

  if (confile.empty()) {
    cerr << colorizer.error(
        "At least one non-option argument is required: the con file\n");
    cerr << colorizer.warning(parser.formatUsage(progname)) << '\n';
    exit(EXIT_FAILURE);
  }

  if (sflag && mflag) {
    cerr << colorizer.error(
        "Cannot specify both minimization and single point\n");
    exit(EXIT_FAILURE);
  }

  if (!pflag && (sflag || mflag)) {
    cerr << colorizer.error("Must specify a potential\n");
    exit(EXIT_FAILURE);
  }

#ifdef WITH_SERVE_MODE
  // Load config file if provided (for potential-specific parameters
  // like model_path, device, length_unit, etc.)
  if (config_path.has_value()) {
    std::ifstream config_file(config_path.value());
    if (!config_file.is_open()) {
      cerr << colorizer.error("Cannot open config file: ")
           << config_path.value() << '\n';
      exit(EXIT_FAILURE);
    }
    params.load(config_path.value());
  }

  // Handle --serve mode (does not require a con file)
  if (serve_spec.has_value()) {
    auto endpoints = parseServeSpec(serve_spec.value());
    if (endpoints.empty()) {
      cerr << colorizer.error("No valid serve endpoints in spec: ")
           << serve_spec.value() << '\n';
      exit(EXIT_FAILURE);
    }
    serveMultiple(endpoints, params);
    exit(EXIT_SUCCESS);
  }

  // Handle -p with serve flags (single potential serve mode)
  if (pflag && !sflag && !mflag && !cflag &&
      (serve_port.has_value() || replicas.has_value() || gateway.has_value())) {
    for (auto &ch : potential) {
      ch = tolower(ch);
    }
    params.potential_options.potential =
        magic_enum::enum_cast<PotType>(potential, magic_enum::case_insensitive)
            .value_or(PotType::UNKNOWN);
    auto host = serve_host.value_or("localhost");
    auto port = serve_port.value_or(12345);
    auto reps = replicas.value_or(1);
    bool gw = gateway.value_or(false);

    if (gw) {
      serveGateway(params, host, port, reps);
    } else if (reps > 1) {
      serveReplicated(params, host, port, reps);
    } else {
      serveMode(params, host, port);
    }
    exit(EXIT_SUCCESS);
  }

  // Config-driven serve (no -p or --serve, just --config with [Serve])
  if (!pflag && !sflag && !mflag && !cflag && config_path.has_value() &&
      !serve_spec.has_value() &&
      (!params.serve_options.endpoints.empty() ||
       params.serve_options.gateway_port > 0 ||
       params.serve_options.replicas > 1)) {
    serveFromConfig(params);
    exit(EXIT_SUCCESS);
  }
#endif

  if (!cflag) {
    for (auto &ch : potential) {
      ch = tolower(ch);
    }
  }

  if (!cflag) {
    params.potential_options.potential =
        magic_enum::enum_cast<PotType>(potential, magic_enum::case_insensitive)
            .value_or(PotType::UNKNOWN);
  }

  if (!sflag) {
    params.optimizer_options.method =
        magic_enum::enum_cast<OptType>(optimizer, magic_enum::case_insensitive)
            .value_or(OptType::CG);
    params.optimizer_options.converged_force = optConvergedForce;
  }

  auto pot = helper_functions::makePotential(params);
  auto matter = std::make_unique<Matter>(pot, params);
  auto matter2 = std::make_unique<Matter>(pot, params);
  matter->con2matter(confile);

  if (confileout.empty() && confile.empty()) {
    // confileout was not provided as second positional arg
    // It's empty, which is fine for single point/minimize
  }

  if (sflag) {
    singlePoint(std::move(matter));
  } else if (mflag) {
    minimize(std::move(matter), confileout);
  } else if (cflag) {
    params.structure_comparison_options.check_rotation = true;
    if (matter->compare(*matter2, true)) {
      cout << "Structures match\n";
    } else {
      cout << colorizer.error("Structures do not match\n");
    }
  }
}
